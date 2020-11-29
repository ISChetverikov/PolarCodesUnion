
#include <algorithm>
#include <sstream>
#include <iostream>
#include "../include/ScListFlipOracleStatDecoder.h"

#define FROZEN_VALUE 0

ScListFlipOracleStatDecoder::ScListFlipOracleStatDecoder(PolarCode * codePtr, int L) : ScListDecoder(codePtr, L) {
	_level = 2;
}

void ScListFlipOracleStatDecoder::ClearStatistic() {
	_singleOracleFlipsStat.clear();
	_doubleOracleFlipsStat.clear();
	_tripleOracleFlipsStat.clear();

	_llrWhenFlip.clear();
	_count = 0;
	_countErrorneous = 0;
	_singleSuccessfulFlips = 0;
	_doubleSuccessfulFlips = 0;
	_tripleSuccessfulFlips = 0;
}

std::string ScListFlipOracleStatDecoder::GetStatistic() {
	std::stringstream ss;
	
	bool saveStat = false;
	bool saveLlrs = false;

	if (saveStat) {
		ss << "Single Flip:\n";
		for (std::pair<int, int> single : _singleOracleFlipsStat)
		{
			ss << "(" << single.first << "): " << single.second << "\n";
		}

		ss << "Double Flip:\n";
		for (std::pair<std::tuple<int, int>, int> double_ : _doubleOracleFlipsStat)
		{
			ss << "(" << std::get<0>(double_.first) << ", " << std::get<1>(double_.first) << "): " << double_.second << "\n";
		}

		ss << "Triple Flip:\n";
		for (std::pair<std::tuple<int, int, int>, int> triple : _tripleOracleFlipsStat)
		{
			ss << "(" << std::get<0>(triple.first) << ", " << std::get<1>(triple.first) << ", " << std::get<2>(triple.first) << "): " << triple.second << "\n";
		}
	}
	
	if (saveLlrs) {
		ss << "Lllrs when flip occurs:\n";
		for (size_t i = 0; i < _llrWhenFlip.size(); i++)
		{
			ss << "(" << std::get<0>(_llrWhenFlip[i]) << ") {";

			for (size_t j = 0; j < _L; j++)
			{
				ss << std::get<1>(_llrWhenFlip[i])[j] << " ";
			}

			ss << "} [";

			for (size_t j = 0; j < 2 * _L; j++)
			{
				if (j == _L)
					ss << "| ";
				ss << std::get<2>(_llrWhenFlip[i])[j] << " ";
			}
			ss << "]\n";
		}
	}
	
	ss << "Count:" << _count << "\n";
	ss << "List Error Count: " << _countErrorneous << "\n";
	ss << "Single: " << _singleSuccessfulFlips << "\n";
	ss << "Double: " << _doubleSuccessfulFlips << "\n";
	ss << "Triple: " << _tripleSuccessfulFlips << "\n";

	return ss.str();
}

void ScListFlipOracleStatDecoder::SaveLlrs(int i_all) {
	std::vector<double> llrs(_L, 0);
	size_t m = _codePtr->m();

	for (size_t j = 0; j < _L; j++)
	{
		llrs[j] = _beliefTrees[j][m][i_all];
	}

	std::vector<double> metrics(2 * _L, 0);
	for (size_t j = 0; j < _L; j++)
	{
		metrics[j] = _metrics[j] + StepMetric(_beliefTrees[j][m][i_all], 0);
		metrics[j + _L] = _metrics[j] + StepMetric(_beliefTrees[j][m][i_all], 1);
	}

	auto tuple = std::make_tuple(i_all, llrs, metrics);
	_llrWhenFlip.push_back(tuple);
}

std::vector<int> ScListFlipOracleStatDecoder::DecodeFlipOracleListInternal(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();

	//inLlr = { 0.94217, 0.918489, 0.891565, 0.0601442, 0.661353, 0.627861, 0.0147878, 0.310293, 0.157784, 0.65148, 0.999905, 0.802843, 0.983206, 0.0200608, 0.0797412, 0.266692 };

	// for saving flips
	std::vector<int> flips;

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

			size_t j = 0;
			size_t i = 0;
			while (i < i_all)
			{
				if (_codeword[i] != _candidates[j][i]) {
					j++;
					i = 0;
				}
				i++;
			}

			int rightBit = _codeword[i_all];
			if (rightBit == 1 && _areTakenOne[j] == 0 || rightBit == 0 && _areTakenZero[j] == 0) {

				if (flips.size() >= _level)
					break;

				for (size_t i = 0; i < _L; i++) {
					_areTakenOne[i] = !_areTakenOne[i];
					_areTakenZero[i] = !_areTakenZero[i];
				}
				flips.push_back((int)i_all);
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

	return flips;
}

std::vector<int>  ScListFlipOracleStatDecoder::Decode(std::vector<double> beliefs) {
	_count++;

	std::vector<int> flipPositions = DecodeFlipOracleListInternal(beliefs);

	size_t k = flipPositions.size();

	if (k > 0)
		_countErrorneous++;

	if (k == 1) {
		_singleSuccessfulFlips++;
		_singleOracleFlipsStat[flipPositions[0]]++;
	}
		

	if (k == 2) {
		_doubleSuccessfulFlips++;
		_doubleOracleFlipsStat[std::make_tuple(flipPositions[0], flipPositions[1])]++;
	}
		

	if (k == 3) {
		_tripleSuccessfulFlips++;
		_tripleOracleFlipsStat[std::make_tuple(flipPositions[0], flipPositions[1], flipPositions[2])]++;
	}
		

	return TakeListResult();
}