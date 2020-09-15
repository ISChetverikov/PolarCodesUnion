#pragma once

#include "ScCrcAidedDecoder.h"

struct CriticalSetNode {
public:
	std::vector<CriticalSetNode *> Children;
	std::vector<int> Path;
	int Bit;
};

class ScFlipProgDecoder : public ScCrcAidedDecoder {

	
private:
	
	CriticalSetNode * _criticalSetTree;
	std::vector<double> _subchannelsMeansGa;

	double _gammaLeft;
	double _gammaRight;
	std::vector<double> _omegaArr;
	int _levelMax;

	std::vector<int> GetCriticalSet(std::vector<int> mask, int position);
	std::vector<CriticalSetNode *> SortCriticalNodes(std::vector<CriticalSetNode *> criticalSet, std::vector<double> llrs);
	CriticalSetNode * GetCriticalSetTree(std::vector<int> mask, int levelMax);
	bool NoChild(CriticalSetNode * node, std::vector<double> inLlr);
	bool NotSelect(CriticalSetNode * node, std::vector<double> inLlr);

public:
	ScFlipProgDecoder(PolarCode * codePtr,
		int level,
		double gammaLeft,
		double gammaRight,
		std::vector<double> omegaArr);

	std::vector<int> Decode(std::vector<double> llr) override;
	void SetSigma(double sigma) override;

	~ScFlipProgDecoder() {};
};