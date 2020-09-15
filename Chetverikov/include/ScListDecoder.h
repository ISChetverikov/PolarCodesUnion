#pragma once

#include "ScCrcAidedDecoder.h"

class ScListDecoder : public ScCrcAidedDecoder {
protected:
	int _L;
	std::vector<std::vector<std::vector<double>>> _beliefTrees;
	std::vector<std::vector<std::vector<int>>> _uhatTrees;
	std::vector<std::vector<int>> _candidates;
	std::vector<double> _metrics;

	// optimization allocation
	std::vector<bool> _areTakenZero;
	std::vector<bool> _areTakenOne;

	void PassDownList(size_t iter);
	void PassUpList(size_t iter);
	void DecodeListInternal(std::vector<double> inLlr);
	void FillListMask(size_t iter);
	double StepMetric(double belief, int decision);
	std::vector<int> TakeListResult();

public:
	ScListDecoder(PolarCode * code, int L);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListDecoder() {};
};