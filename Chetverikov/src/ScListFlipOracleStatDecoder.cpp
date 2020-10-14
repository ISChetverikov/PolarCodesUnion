
#include <sstream>
#include "../include/ScListFlipOracleStatDecoder.h"

#define FROZEN_VALUE 0

ScListFlipOracleStatDecoder::ScListFlipOracleStatDecoder(PolarCode * codePtr, int L) : ScListFlipStatDecoder(codePtr, L) {
}

int ScListFlipOracleStatDecoder::GetFirstErrorPosition(std::vector<int> codeword1, std::vector<int> codeword2) {
	int i = 0;

	for (; i < codeword1.size(); i++)
	{
		if (codeword1[i] != codeword2[i])
			break;
	}

	return i;
}

void ScListFlipOracleStatDecoder::ClearStatistic() {
	_singleOracleFlipsStat.clear();
	_doubleOracleFlipsStat.clear();
	_tripleOracleFlipsStat.clear();
}

std::string ScListFlipOracleStatDecoder::GetStatistic() {
	std::stringstream ss;
	// std::sort(_singleFlipStatistic.rbegin(), _singleFlipStatistic.rend());
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

	return ss.str();
}



std::vector<int>  ScListFlipOracleStatDecoder::Decode(std::vector<double> beliefs) {
	std::vector<int> actualCodeword(beliefs.size(), 0);
	std::vector<int> flipPositions;
	bool isError = false;
	
	DecodeListInternal(beliefs);
		
	std::vector<int> result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError)
		return result;

	int firstErrorInd = GetFirstErrorPosition(actualCodeword, _codeword);
	flipPositions.clear();
	flipPositions.push_back(firstErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);
	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		if (_singleOracleFlipsStat.find(firstErrorInd) == _singleOracleFlipsStat.end())
			_singleOracleFlipsStat[firstErrorInd] = 1;
		else
			_singleOracleFlipsStat[firstErrorInd]++;

		return result;
	}
	
	int secondErrorInd = GetFirstErrorPosition(actualCodeword, _codeword);
	flipPositions.push_back(secondErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);
	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		auto tuple = std::make_tuple(firstErrorInd, secondErrorInd);

		if (_doubleOracleFlipsStat.find(tuple) == _doubleOracleFlipsStat.end())
			_doubleOracleFlipsStat[tuple] = 1;
		else
			_doubleOracleFlipsStat[tuple]++;

		return result;
	}
	
	int thirdErrorInd = GetFirstErrorPosition(actualCodeword, _codeword);
	flipPositions.push_back(thirdErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);
	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		auto tuple = std::make_tuple(firstErrorInd, secondErrorInd, thirdErrorInd);

		if (_tripleOracleFlipsStat.find(tuple) == _tripleOracleFlipsStat.end())
			_tripleOracleFlipsStat[tuple] = 1;
		else
			_tripleOracleFlipsStat[tuple]++;

		return result;
	}

	return result;
}