#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "PolarCode.h"
#include "Domain.h"

class BaseDecoder {

protected:
	PolarCode * _codePtr;
	double _sigma;
public:
	BaseDecoder(PolarCode * codePtr);
	
	virtual domain GetDomain() = 0;
	void SetSigma(double sigma);
	double GetSigma();
	virtual std::vector<int> Decode(std::vector<double> llr) = 0;
	virtual ~BaseDecoder() {};
};