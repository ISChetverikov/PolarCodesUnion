#pragma once

#include <string>
#include <math.h>
#include "Encoder.h"
#include "PolarCode.h"
#include "BaseDecoder.h"
#include "SimulationParameters.h"
#include "BaseChannel.h"

class BaseSimulator {
protected:
	PolarCode * _codePtr;
	Encoder * _encoderPtr;
	BaseChannel * _channelPtr;
	BaseDecoder * _decoderPtr;
	bool _isSigmaDependOnR;
	std::string _additionalInfoFilename;

	double GetSigma(double snr, double R);
	double GetEbN0(double snr, size_t m, size_t n);
public:
	BaseSimulator(PolarCode * codePtr, Encoder * encoder, BaseChannel * channelPtr, BaseDecoder * decoder, bool isEbnoSimulation);
	virtual ~BaseSimulator() {};
	virtual SimulationIterationResults Run(double snr) = 0;
};