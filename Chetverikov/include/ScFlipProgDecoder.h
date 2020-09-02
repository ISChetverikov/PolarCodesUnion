#pragma once

#include "BaseDecoder.h"
#include "Domain.h"
#include "CRC.h"

struct CriticalSetNode {
public:
	std::vector<CriticalSetNode *> Children;
	std::vector<int> Path;
	int Bit;
};

class ScFlipProgDecoder : public BaseDecoder {

	
private:
	std::vector<std::vector<double>> _beliefTree;
	std::vector<std::vector<int>> _uhatTree;
	size_t _treeHeight;
	std::vector<int> _maskWithCrc;
	std::vector<int> _x;
	CRC * _crcPtr;
	CriticalSetNode * _criticalSetTree;
	std::vector<double> _subchannelsMeansGa;

	double _gammaLeft;
	double _gammaRight;
	std::vector<double> _omegaArr;
	int _levelMax;

	// optimization allocation
	std::vector<int> _binaryIter;

	double f(double llr1, double llr2);
	double g(double llr1, double llr2, int u1);
	int L(double llr1);
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

	void DecodeFromTo(int position, int endPosition);
	bool IsCrcPassed(std::vector<int> codeword);

	std::vector<int> GetCriticalSet(std::vector<int> mask, int position);
	std::vector<CriticalSetNode *> SortCriticalNodes(std::vector<CriticalSetNode *> criticalSet, std::vector<double> llrs);
	CriticalSetNode * GetCriticalSetTree(std::vector<int> mask, int levelMax);
	bool NoChild(CriticalSetNode * node, std::vector<double> inLlr);
	bool NotSelect(CriticalSetNode * node, std::vector<double> inLlr);

public:
	ScFlipProgDecoder(PolarCode * code, int level, double gammaLeft, double gammaRight, std::vector<double> omegaArr);

	std::vector<int> Decode(std::vector<double> llr) override;
	void SetSigma(double sigma) override;
	domain GetDomain() override;

	~ScFlipProgDecoder() {};
};