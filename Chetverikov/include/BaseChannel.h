#pragma once

#include <vector>

class BaseChannel {
protected:
	double _snr = 1.0;

public:
	BaseChannel();
	virtual ~BaseChannel() {};
	virtual void SetSnr(double sigma);
	virtual std::vector<double> Pass(std::vector<int> codeword) = 0;
};