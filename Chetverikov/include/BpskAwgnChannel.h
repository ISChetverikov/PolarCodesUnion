#pragma once

#include <random>
#include "BaseChannel.h"

class BpskAwgnChannel : public BaseChannel {
protected:
	std::random_device _randomDevice;
	std::normal_distribution<double> _normal_dist;
	double _sigma;
public:
	BpskAwgnChannel();
	~BpskAwgnChannel() {};
	std::vector<double> Pass(std::vector<int> codeword) override;
	void SetSnr(double sigma) override;
};