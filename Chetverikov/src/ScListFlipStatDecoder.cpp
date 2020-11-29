#include <sstream>
#include <algorithm>
#include <iostream>
#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScListFlipStatDecoder.h"
#include "../include/FlipStatisticNew256_128.h"

#define FROZEN_VALUE 0

ScListFlipStatDecoder::ScListFlipStatDecoder(PolarCode * codePtr, int L, double omega) : ScListDecoder(codePtr, L) {
	//_unfrozenPolarSeqWithCrc = _codePtr->UnfrozenPolarSequenceWithCrc();
	int singleNumber = 40; // TODO
	int doubleNumber = -1;
	bool isUsedStat = true;
	bool isDoubleFromSingles = false;

	if (isUsedStat) {
		_singleFlips = FlipStatistic::GetSingles();

		if (isDoubleFromSingles) {
			_doubleFlips = std::vector<std::tuple<int, int>>();

			for (size_t i = 0; i < _singleFlips.size(); i++)
			{
				for (size_t j = i + 1; j < _singleFlips.size(); j++)
				{
					_doubleFlips.push_back(std::make_tuple(_singleFlips[i], _singleFlips[j]));
				}
			}
		}
		else
			_doubleFlips = FlipStatistic::GetPairs();
	}
	else {
		std::vector<int> unfrozenBitsWithCrc = codePtr->UnfrozenPolarSequenceWithCrc();

		_singleFlips = unfrozenBitsWithCrc;
		_doubleFlips = std::vector<std::tuple<int, int>>();
		for (size_t i = 0; i < unfrozenBitsWithCrc.size(); i++)
		{
			for (size_t j = i + 1; j < unfrozenBitsWithCrc.size(); j++)
			{
				_doubleFlips.push_back(std::make_tuple(unfrozenBitsWithCrc[i], unfrozenBitsWithCrc[j]));
			}
		}
	}

	singleNumber = _singleFlips.size();
	_singleFlips = std::vector<int>(_singleFlips.begin(), _singleFlips.begin() + singleNumber);
	_omega = omega;
}

bool ScListFlipStatDecoder::IsThresholdPassed(int i_all) {
	return false;

	// turned off
	auto maxIt = std::max_element(_metrics.begin(), _metrics.end());
	int maxInd = (int)std::distance(_metrics.begin(), maxIt);

	size_t m = _codePtr->m();
	return fabs(_beliefTrees[maxInd][m][i_all]) > _omega;
}

void ScListFlipStatDecoder::ClearStatistic() {
	_countThreasholdTries = 0;
	_countThresholdNotPassed = 0;
}

std::string ScListFlipStatDecoder::GetStatistic() {
	std::stringstream ss;

	ss << _countThresholdNotPassed << " from " << _countThreasholdTries << ": " << (double)_countThresholdNotPassed / _countThreasholdTries << "\n";

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
			if (std::find(flipPositions.begin(), flipPositions.end(), i_all) != flipPositions.end()) {
				_countThreasholdTries++;
				
				if (IsThresholdPassed(i_all)) {
					break;
				}
				_countThresholdNotPassed++;

				for (size_t i = 0; i < _L; i++) {
					_areTakenOne[i] = !_areTakenOne[i];
					_areTakenZero[i] = !_areTakenZero[i];
				}
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

std::vector<int> ScListFlipStatDecoder::TakeListStatResult(bool& isError) {
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
	if (j == _L) {
		isError = true;
		return std::vector<int>();
	}

	for (size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = candidate[codewordBits[i]];
	}

	isError = false;
	return result;
}

std::vector<int>  ScListFlipStatDecoder::Decode(std::vector<double> beliefs) {
	DecodeListInternal(beliefs);
	
	bool isError = false;
	std::vector<int> result = TakeListStatResult(isError);

	if (!isError)
		return result;

	std::vector<int> flipPositions;

	for (size_t i = 0; i < _singleFlips.size(); i++)
	{
		flipPositions.clear();
		flipPositions.push_back(_singleFlips[i]);

		DecodeFlipListInternal(beliefs, flipPositions);
		result = TakeListStatResult(isError);

		if (!isError) {
			return result;
		}
	}

	return result;
}