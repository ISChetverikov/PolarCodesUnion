#pragma once

#include "ScListDecoder.h"

class ScListFlipStatDecoder : public ScListDecoder {
private:
	std::vector<int> _unfrozenPolarSeqWithCrc;
	std::vector<int> _singleFlipStatistic;
	std::vector<std::vector<int>> _doubleFlipStatistic;

	void DecodeFlipListInternal(std::vector<double> inLlr, std::vector<int> flipPositions);
	std::vector<int> TakeListStatResult(bool & isOriginalCodeword);
public:
	ScListFlipStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipStatDecoder() {};
};