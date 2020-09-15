#pragma once

#include "CRC.h"
#include "ScDecoder.h"

// logically it is abstract class
class ScCrcAidedDecoder : public ScDecoder {

protected:
	
	CRC * _crcPtr;

	void DecodeFrom(int position);
	void DecodeFromTo(int position, int endPosition);
	bool IsCrcPassed(std::vector<int> codeword);

public:
	ScCrcAidedDecoder(PolarCode * code);
	
	~ScCrcAidedDecoder() {};
};