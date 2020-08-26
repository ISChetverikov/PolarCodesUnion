#pragma once

#include "PolarCode.h"

class Encoder {

protected:
	PolarCode * _codePtr;

public:
	Encoder();
	Encoder(PolarCode * codePtr);

	std::vector<int> Encode(std::vector<int> x);
	~Encoder();
};