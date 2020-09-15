#pragma once



#include "ScDecoder.h"

class ScFanoDecoder : public ScDecoder {

private:
	double _T;
	double _delta;

	std::vector<double> _p; // channel error probabilities 

public:
	ScFanoDecoder(PolarCode * code, double T, double delta);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFanoDecoder() {};
};