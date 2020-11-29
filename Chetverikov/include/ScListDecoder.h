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
	// At each step: _areTakenZero[i] = 1 if i-th path with added zero is taken
	std::vector<bool> _areTakenZero;
	std::vector<bool> _areTakenOne;

	void PassDownList(size_t iter);
	void PassUpList(size_t iter);
	void DecodeListInternal(std::vector<double> inLlr);
	// fill _areTakenZero and _areTakeOne
	void FillListMask(size_t iter);
	double StepMetric(double belief, int decision);
	std::vector<int> TakeListResult();

public:
	ScListDecoder(PolarCode * code, int L);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListDecoder() {};
};