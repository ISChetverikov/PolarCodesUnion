#pragma once

#include "ScListDecoder.h"

class ScListFlipStatDecoder : public ScListDecoder {
private:
	std::vector<int> _singleFlips;
	std::vector<std::tuple<int, int>> _doubleFlips;

	std::vector<int> _singleFlipStatistic;
	std::vector<int> _doubleFlipStatistic;

protected:
	void DecodeFlipListInternal(std::vector<double> inLlr, std::vector<int> flipPositions);
	std::vector<int> TakeListStatResult(std::vector<int> & actualCodeword);
public:
	ScListFlipStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipStatDecoder() {};
};