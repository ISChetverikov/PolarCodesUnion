#pragma once

#include "ScCrcAidedDecoder.h"

class ScFlipDecoder : public ScCrcAidedDecoder {
private:
	int _T;
	std::vector<int> GetSmallestBeliefsIndices(std::vector<double> llrs, int count);

public:
	ScFlipDecoder(PolarCode * code, int T);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFlipDecoder() {};
};