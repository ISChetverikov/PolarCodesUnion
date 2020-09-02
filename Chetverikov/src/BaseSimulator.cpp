#include "../include/BaseSimulator.h"

BaseSimulator::BaseSimulator(PolarCode * codePtr, Encoder * encoderPtr, BaseDecoder * decoderPtr, bool isSigmaDependOnR) {
	_decoderPtr = decoderPtr;
	_encoderPtr = encoderPtr;
	_codePtr = codePtr;
	_isSigmaDependOnR = isSigmaDependOnR;
}

double BaseSimulator::GetSigma(double snr, double R) {
	double result = sqrt(pow(10, -snr / 10) / 2);
	if (_isSigmaDependOnR)
		result /= R;

	return result;
}

double BaseSimulator::GetEbN0(double snr, size_t m, size_t n) {
	return snr - 10 * log10((double)m / n);
}