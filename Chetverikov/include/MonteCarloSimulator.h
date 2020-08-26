#pragma once

#include <string>
#include <chrono>
#include "Encoder.h"
#include "BaseDecoder.h"
#include "BaseSimulator.h"

class MonteCarloSimulator : public BaseSimulator {
protected:
	int _maxTestsCount = 0;
	int _maxRejectionsCount = 0;
	

public:
	MonteCarloSimulator(int maxTests, int rejectionsCount, PolarCode * codePtr, Encoder * encoderPtr, BaseDecoder * decoderPtr);
	~MonteCarloSimulator() {};
	SimulationIterationResults Run(double snr) override;
};