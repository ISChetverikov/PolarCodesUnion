#pragma once

#include "ScCrcAidedDecoder.h"

class ScFlipFanoDecoder : public ScCrcAidedDecoder {

private:
	double _T;
	double _delta;
	double _L;

	std::vector<double> _p; // channel error probabilities 

	// inner state
	std::vector<double> _beta; // only for unfrozen bits
	std::vector<double> _metrics; // for all bits
	std::vector<double> _alternativeBeta; // unfrozenBits alternative metric
	std::vector<double> _betaDifference; // unfrozenBits alternative metric
	std::vector<bool> _gamma;
	std::vector<int> _A; // info set

	void UpdateT(double & T, double & tau);
	void BackwardMove(double & T, bool & B, int & j, int rootIndex);
	void DecodeFrom(int rootIndex);
	void Flip(int i);

public:
	ScFlipFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, double L);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFlipFanoDecoder() {};
};