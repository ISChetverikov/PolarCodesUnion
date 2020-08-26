#include "../include/BaseSimulator.h"

BaseSimulator::BaseSimulator(PolarCode * codePtr, Encoder * encoderPtr, BaseDecoder * decoderPtr) {
	_decoderPtr = decoderPtr;
	_encoderPtr = encoderPtr;
	_codePtr = codePtr;
}

double BaseSimulator::GetSigma(double snr, double R) {
	return sqrt(pow(10, -snr / 10) / 2 / R);
}

double BaseSimulator::GetEbN0(double snr, size_t m, size_t n) {
	return snr - 10 * log10((double)m / n);
}