#include <sstream>
#include <algorithm>

#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScListFlipStatDecoder.h"

#define FROZEN_VALUE 0

ScListFlipStatDecoder::ScListFlipStatDecoder(PolarCode * codePtr, int L) : ScListDecoder(codePtr, L) {
	_unfrozenPolarSeqWithCrc = _codePtr->UnfrozenPolarSequenceWithCrc();
	ClearStatistic();
}

void ScListFlipStatDecoder::ClearStatistic() {
	size_t kExt = _codePtr->kExt();
	_singleFlipStatistic = std::vector<int>(kExt, 0);
	_doubleFlipStatistic = std::vector<std::vector<int>>(kExt, _singleFlipStatistic);
}

std::string ScListFlipStatDecoder::GetStatistic() {
	std::stringstream ss;
	// std::sort(_singleFlipStatistic.rbegin(), _singleFlipStatistic.rend());
	ss << "Single Flip:\n";
	for (size_t i = 0; i < _codePtr->kExt(); i++)
	{
		// if (_singleFlipStatistic[i])
			ss << "(" << _unfrozenPolarSeqWithCrc[i] << "): " << _singleFlipStatistic[i] << "\n";
	}

	ss << "Double Flip:\n";
	for (size_t i = 0; i < _codePtr->kExt(); i++)
	{
		for (size_t j = 0; j < _codePtr->kExt(); j++)
		{
			std::sort(_doubleFlipStatistic[i].rbegin(), _doubleFlipStatistic[i].rend());
			// if (_doubleFlipStatistic[i][j])
				ss << "(" << _unfrozenPolarSeqWithCrc[i] << ", " << _unfrozenPolarSeqWithCrc[j] << "): " << _doubleFlipStatistic[i][j] << "\n";
			
		}
	}

	return ss.str();
}

void ScListFlipStatDecoder::DecodeFlipListInternal(std::vector<double> inLlr, std::vector<int> flipPositions) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();

	// Fill each tree in the forrest with input llrs
	for (size_t j = 0; j < _L; j++) {
		for (size_t i = 0; i < n; i++)
			_beliefTrees[j][0][i] = inLlr[i];

		_metrics[j] = 0;
	}

	int logL = (int)FirstBitPos(_L) - 1;
	// first log(_L) bits
	size_t i_all = 0; // number of bit (j - number of condidate in list)
	size_t i_unfrozen = 0;
	while (i_unfrozen < logL)
	{
		PassDownList(i_all);

		if (_maskWithCrc[i_all]) {
			int value = 0;
			for (size_t j = 0; j < _L; j++) {
				_candidates[j][i_all] = value;
				if ((j + 1) % (1 << i_unfrozen) == 0) // all paths at the logL first steps 
					value = !value;
			}
			i_unfrozen++;
		}
		else {
			for (size_t j = 0; j < _L; j++)
				_candidates[j][i_all] = FROZEN_VALUE;
		}

		for (size_t j = 0; j < _L; j++) {
			_uhatTrees[j][m][i_all] = _candidates[j][i_all];
			_metrics[j] += StepMetric(_beliefTrees[j][m][i_all], _candidates[j][i_all]);
		}

		PassUpList(i_all);

		i_all++;
	}
	auto c = _codeword;
	while (i_all < n)
	{
		PassDownList(i_all);

		if (_maskWithCrc[i_all]) {
			FillListMask(i_all);

			// Flip action
			if (std::find(flipPositions.begin(), flipPositions.end(), i_all) != flipPositions.end())
				for (size_t i = 0; i < _L; i++) {
					_areTakenOne[i] = !_areTakenOne[i];
					_areTakenZero[i] = !_areTakenZero[i];
				}

			for (size_t j = 0, i = 0; i < _L; i++) {
				// add new path
				if (_areTakenOne[j] && _areTakenZero[j]) {
					_candidates[j][i_all] = 1;
					_beliefTrees.push_back(_beliefTrees[j]);
					_metrics.push_back(_metrics[j]);
					_uhatTrees.push_back(_uhatTrees[j]);
					_candidates.push_back(_candidates[j]);

					_candidates[j][i_all] = 0;
					_areTakenOne.push_back(false);
					_areTakenZero.push_back(false);

					j++;
					continue;
				}
				if (_areTakenZero[j]) {
					_candidates[j][i_all] = 0;
					j++;
					continue;
				}
				if (_areTakenOne[j]) {
					_candidates[j][i_all] = 1;
					j++;
					continue;
				}

				// delete path
				_beliefTrees.erase(_beliefTrees.begin() + j);
				_metrics.erase(_metrics.begin() + j);
				_uhatTrees.erase(_uhatTrees.begin() + j);
				_candidates.erase(_candidates.begin() + j);
				_areTakenZero.erase(_areTakenZero.begin() + j);
				_areTakenOne.erase(_areTakenOne.begin() + j);
			}
		}
		else {
			for (size_t j = 0; j < _L; j++)
				_candidates[j][i_all] = FROZEN_VALUE;
		}

		for (size_t j = 0; j < _L; j++) {
			_uhatTrees[j][m][i_all] = _candidates[j][i_all];
			_metrics[j] += StepMetric(_beliefTrees[j][m][i_all], _candidates[j][i_all]);
		}

		PassUpList(i_all);

		i_all++;
	}

	return;
}

std::vector<int> ScListFlipStatDecoder::TakeListStatResult(bool & isError) {
	std::vector<int> result(_codePtr->k(), 0);
	std::vector<int> candidate(_codePtr->N(), 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();

	int maxInd = -1;
	size_t j = 0;
	for (; j < _L; j++)
	{
		auto maxIt = std::max_element(_metrics.begin(), _metrics.end());
		maxInd = (int)std::distance(_metrics.begin(), maxIt);

		if (IsCrcPassed(_candidates[maxInd]))
			break;

		_metrics[maxInd] = -10000000.0;
	}

	if (j < _L) {
		candidate = _candidates[maxInd];	
	}
	isError = (candidate != _codeword);

	for (size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = candidate[codewordBits[i]];
	}

	return result;
}

std::vector<int>  ScListFlipStatDecoder::Decode(std::vector<double> beliefs) {
	DecodeListInternal(beliefs);

	bool isError = false;
	std::vector<int> result = TakeListStatResult(isError);

	if (!isError)
		return result;

	size_t kExt = _codePtr->kExt();
	for (size_t i = 0; i < kExt; i++)
	{
		int bitPosition = _unfrozenPolarSeqWithCrc[i];

		std::vector<int> flipPositions;
		flipPositions.push_back(bitPosition);

		DecodeFlipListInternal(beliefs, flipPositions);
		result = TakeListStatResult(isError);

		if (!isError) {
			_singleFlipStatistic[i]++;
			return result;
		}
	}
	
	for (size_t i = 0; i < kExt; i++)
	{
		for (size_t j = i + 1; j < kExt; j++)
		{
			int bitPosition1 = _unfrozenPolarSeqWithCrc[i];
			int bitPosition2 = _unfrozenPolarSeqWithCrc[j];

			std::vector<int> flipPositions;
			flipPositions.push_back(bitPosition1);
			flipPositions.push_back(bitPosition2);

			DecodeFlipListInternal(beliefs, flipPositions);
			result = TakeListStatResult(isError);

			if (!isError) {
				_doubleFlipStatistic[i][j]++;
				return result;
			}
		}
	}

	return result;

	/*size_t k = _codePtr->k();
	for (size_t i = 0; i < k; i++)
	{
		int bitPosition = _unfrozenBits[i];

		if (_x[bitPosition] == originalCodeword[bitPosition])
			continue;

		_x[bitPosition] = !_x[bitPosition];

		DecodeFrom(bitPosition);

		if (_x == _codeword) {
			_singleFlipStatistic[i]++;
			return TakeResult();
		}

		for (size_t j = i + 1; j < k; j++)
		{
			int bitPosition2 = _unfrozenBits[j];
			if (_x[bitPosition2] == originalCodeword[bitPosition2])
				continue;

			_x[bitPosition2] =! _x[bitPosition2];
			DecodeFrom(bitPosition2);

			if (_x == _codeword) {
				_doubleFlipStatistic[i][j]++;
				return TakeResult();
			}

			_x[bitPosition2] = !_x[bitPosition2];
			DecodeFrom(bitPosition2);
		}

		_x[bitPosition] = !_x[bitPosition];
		DecodeFrom(bitPosition);
	}*/
}