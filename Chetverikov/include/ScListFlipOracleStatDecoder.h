#pragma once

#include <map>
#include "ScListDecoder.h"

class ScListFlipOracleStatDecoder : public ScListDecoder {
private:
	std::map<int, int> _singleOracleFlipsStat;
	std::map<std::tuple<int, int>, int> _doubleOracleFlipsStat;
	std::map<std::tuple<int, int, int>, int> _tripleOracleFlipsStat;

	// want to learn llr when flipping is occured
	// tuple: flip position: int, vector of llr, vector metrics: the first part - 0 step, the second - 1 flip
	std::vector<std::tuple<int, std::vector<double>, std::vector<double>>> _llrWhenFlip;

	int _count = 0;
	int _countErrorneous = 0;
	int _singleSuccessfulFlips = 0;
	int _doubleSuccessfulFlips = 0;
	int _tripleSuccessfulFlips = 0;
	
	int _level = 0;

	void SaveLlrs(int i);
protected:

	// Returns exact flips by comparision with _codeword at each step of list decoding
	std::vector<int> DecodeFlipOracleListInternal(std::vector<double> inLlr);

public:
	ScListFlipOracleStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipOracleStatDecoder() {};
};