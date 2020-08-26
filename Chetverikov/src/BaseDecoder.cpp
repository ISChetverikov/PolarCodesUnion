#include "../include/Exceptions.h"
#include "../include/BaseDecoder.h"
#include "../include/Domain.h"

BaseDecoder::BaseDecoder(PolarCode * codePtr) {
	_codePtr = codePtr;
	_sigma = 0;
}

double BaseDecoder::GetSigma() {
	return _sigma;
}

void BaseDecoder::SetSigma(double sigma) {
	_sigma = sigma;
}