#pragma once

#include "BaseDecoder.h"

class ScDecoder : public BaseDecoder {

protected:
	std::vector<std::vector<double>> _beliefTree;
	std::vector<std::vector<int>> _uhatTree;
	size_t _treeHeight;
	std::vector<int> _maskWithCrc;
	std::vector<int> _x;

	// optimization allocation
	std::vector<int> _binaryIter;

	size_t FirstBitPos(size_t n);

	void FillLeftMessageInTree(std::vector<double>::iterator leftIt,
		std::vector<double>::iterator rightIt,
		std::vector<double>::iterator outIt,
		size_t n);
	void FillRightMessageInTree(std::vector<double>::iterator leftIt,
		std::vector<double>::iterator rightIt,
		std::vector<int>::iterator uhatIt,
		std::vector<double>::iterator outIt,
		size_t n);

	void PassDown(size_t iter);
	void PassUp(size_t iter);
	void DecodeInternal(std::vector<double> inLlr);
	std::vector<int> TakeResult();

public:
	ScDecoder(PolarCode * code);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScDecoder() {};
};