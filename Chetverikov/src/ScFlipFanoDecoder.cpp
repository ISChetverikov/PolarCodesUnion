
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScFlipFanoDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScFlipFanoDecoder::ScFlipFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, double L) : ScCrcAidedDecoder(code) {
	size_t n = _codePtr->N();
	size_t k = _codePtr->kExt();

	_T = T;
	_delta = delta;
	_L = L;

	double approximationSigma = sqrt(pow(10, -approximationSnr / 10.0));
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}

	// Inners states
	_beta = std::vector<double>(k, 0.0); // only for unfrozen bits
	_alternativeBeta = std::vector<double>(k, 0.0);
	_metrics = std::vector<double>(n, 0); // for all bits
	_betaDifference = std::vector<double>(k, 0.0);
	_gamma = std::vector<bool>(k, 0);
	_A = _codePtr->UnfrozenBitsWithCrc(); // info set
}

void ScFlipFanoDecoder::UpdateT(double & T, double & tau) {
	while (T + _delta < tau)
		T += _delta;
}

void ScFlipFanoDecoder::BackwardMove(double & T, bool & B, int & j, int rootIndex) {

	while (true) {
		double mu = 0;

		if (j <= rootIndex)
			mu = -1000;

		if (j >= 1) // or the rootIndex, mmm?
			mu = _beta[j - 1];

		if (mu >= T) {
			j--;
			if (!_gamma[j + 1])
			{
				B = true;
				return;
			}
		}
		else {
			T -= _delta;
			B = false;
			return;
		}

	}
}

// rootIndex - index before the root (e.g. -1 for 0)
void ScFlipFanoDecoder::DecodeFrom(int rootIndex) {
	size_t n = _codePtr->N();
	size_t m = _codePtr->m();
	size_t k = _codePtr->kExt();
	
	int j = rootIndex;
	int i = 0;
	if (j >= 0)
		i = _A[j] + 1;

	bool B = false;
	double T = _T;

	while (i < n)
	{
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][i];
		double p1 = _beliefTree[m][i];

		if (_maskWithCrc[i]) {
			double previous = (i == 0) ? 0 : _metrics[i - 1];
			double m0 = previous + log(p0 / (1 - _p[i]));
			double m1 = previous + log(p1 / (1 - _p[i]));

			double max = (m1 > m0) ? m1 : m0;
			int argmax = (m1 > m0) ? 1 : 0;

			double min = (m1 > m0) ? m0 : m1;
			if (min < -10000.0)
				min = -10000.0;

			int argmin = (m1 > m0) ? 0 : 1;

			if (max > T) {
				if (!B) {
					_x[i] = argmax;
					_uhatTree[m][i] = _x[i];
					PassUp(i);

					_metrics[i] = max;
					_beta[j + 1] = max;
					_alternativeBeta[j + 1] = min;
					_gamma[j + 1] = false;

					double mu = 0;
					if (j != rootIndex)
						mu = _beta[j];

					if (mu < T + _delta)
						UpdateT(T, _beta[j + 1]);
					i++;
					j++;
				}
				else {
					if (min > T) {
						_x[i] = argmin;
						_uhatTree[m][i] = _x[i];
						PassUp(i);

						_metrics[i] = min;
						_beta[j + 1] = min;
						_gamma[j + 1] = true;
						B = false;

						i++;
						j++;
					}
					else {
						if (j == rootIndex) {
							T = T - _delta;
							B = false;
						}
						else {
							BackwardMove(T, B, j, rootIndex);
							i = _A[j + 1];
						}
					}
				}
			}
			else {
				if (j == rootIndex)
					T = T - _delta;

				else {
					BackwardMove(T, B, j, rootIndex);
					i = _A[j + 1];
				}

			}
		}
		else {
			_x[i] = FROZEN_VALUE;
			_uhatTree[m][i] = _x[i];
			PassUp(i);

			double currentMetric = log(p0) - log(1 - _p[i]);
			// cumulative
			_metrics[i] = currentMetric;
			if (i != 0)
				_metrics[i] += _metrics[i - 1];

			i++;
		}
	}
}

void ScFlipFanoDecoder::Flip(int j) {
	size_t m = _codePtr->m();

	int i = _A[j];
	_x[i] = !_x[i];
	_metrics[i] = _alternativeBeta[j];
	_beta[j] = _alternativeBeta[j];
	_uhatTree[m][i] = _x[i];
	PassUp(i);
}

std::vector<int> ScFlipFanoDecoder::Decode(std::vector<double> beliefs) {
	size_t n = beliefs.size();
	size_t k = _codePtr->kExt();
	size_t m = _codePtr->m();

	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = beliefs[i];
	}

	int rootIndex = -1; // j started position
	std::vector<int> flipPositions; // array of j which are needed to be flipped
	std::vector<double> candidatesMetric;

	DecodeFrom(rootIndex);

	if (IsCrcPassed(_x))
		return TakeResult();

	for (size_t i = 0; i < k; i++)
	{
		_betaDifference[i] = fabs(_beliefTree[m][_A[i]] - 0.5);
	}
	
	for (size_t j = 0; j < _L - 1; j++)
	{
		double minBelief = 100000.0;
		int minInd = -1;

		for (size_t i = 0; i < k; i++)
		{
			if (_betaDifference[i] < minBelief) {
				minBelief = _betaDifference[i];
				minInd = i;
			}
		}

		if (minInd == -1)
			break;

		_betaDifference[minInd] = 100000.0;

		flipPositions.push_back(minInd);
	}

	for (size_t l = 0; l < flipPositions.size(); l++)
	{
		Flip(flipPositions[l]);
		DecodeFrom(flipPositions[l]);
		
		if (IsCrcPassed(_x))
			return TakeResult();
		
		Flip(flipPositions[l]);
		DecodeFrom(flipPositions[l]);
	}
	   

	return TakeResult();
}