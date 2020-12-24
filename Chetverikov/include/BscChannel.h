#pragma once

#include "BaseChannel.h"
#include <vector>

class BscChannel : public BaseChannel {
protected:
	std::random_device _randomDevice;
	std::bernoulli_distribution _bernoulli_dist;
	double _coderate;
	double _fixedLlr;

public:
	BscChannel(double coderate);
	~BscChannel() {};
	std::vector<double> Pass(std::vector<int> codeword) override;
	void SetSnr(double snr) override;
};