#pragma once

#include "PolarCode.h"
#include "CRC.h"

class Encoder {

protected:
	PolarCode * _codePtr;
	CRC * _crcPtr;
public:
	Encoder();
	Encoder(PolarCode * codePtr);

	std::vector<int> Encode(std::vector<int> x);
	~Encoder();
};