#pragma once



#include "BaseDecoder.h"

class ScFanoDecoder : public BaseDecoder {

private:
	std::vector<std::vector<double>> _beliefTree;
	std::vector<std::vector<int>> _uhatTree;
	size_t _treeHeight;
	std::vector<int> _mask;
	std::vector<int> _x;
	double _T;
	double _delta;

	std::vector<double> _p; // channel error probabilities 

	double f(double p1, double p2);
	double g(double p1, double p2);
	int L(double p1);
	size_t log2(size_t n);

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

public:
	ScFanoDecoder(PolarCode * code, double T, double delta);
	
	domain GetDomain() override;
	std::vector<int> Decode(std::vector<double> llr) override;

	~ScFanoDecoder() {};
};