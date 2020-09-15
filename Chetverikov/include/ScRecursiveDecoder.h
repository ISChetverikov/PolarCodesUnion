#pragma once

#include "BaseDecoder.h"
#include "Domain.h"

class ScRecursiveDecoder : public BaseDecoder {

protected:
	std::vector<double> GetLeftMessage(std::vector<double>::iterator leftIt, std::vector<double>::iterator rightIt, size_t n);

	std::vector<double> GetRightMessage(std::vector<double>::iterator leftIt, std::vector<double>::iterator rightIt, std::vector<int> uHat, size_t n);

	std::vector<int> DecodeRecursive(
		std::vector<double>::iterator inLlrIt,
		std::vector<double>::iterator outLlrIt,
		std::vector<int>::iterator maskIt,
		std::vector<int>::iterator xIt,
		size_t n
	);

public:
	ScRecursiveDecoder(PolarCode * codePtr);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScRecursiveDecoder() {};
};