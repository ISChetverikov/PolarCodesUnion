#pragma once

#include <string>
#include <math.h>
#include "Encoder.h"
#include "PolarCode.h"
#include "BaseDecoder.h"
#include "SimulationParameters.h"

class BaseSimulator {
protected:
	PolarCode * _codePtr;
	Encoder * _encoderPtr;
	BaseDecoder * _decoderPtr;
	bool _isSigmaDependOnR;

	double GetSigma(double snr, double R);
	double GetEbN0(double snr, size_t m, size_t n);
public:
	BaseSimulator(PolarCode * codePtr, Encoder * encoder, BaseDecoder * decoder, bool isSigmaDependOnR);
	virtual ~BaseSimulator() {};
	virtual SimulationIterationResults Run(double snr) = 0;
};