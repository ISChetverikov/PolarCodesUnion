#pragma once

#include "ScListDecoder.h"

class ScListFlipStatDecoder : public ScListDecoder {
private:
	std::vector<int> _singleFlips;
	std::vector<std::tuple<int, int>> _doubleFlips;

	double _omega;
	int _countThreasholdTries = 0 ;
	int _countThresholdNotPassed = 0;
protected:
	void DecodeFlipListInternal(std::vector<double> inLlr, std::vector<int> flipPositions);
	std::vector<int> TakeListStatResult(bool& isError);

	bool IsThresholdPassed(int iterationNumber);
public:
	ScListFlipStatDecoder(PolarCode * code, int L, double omega);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipStatDecoder() {};
};