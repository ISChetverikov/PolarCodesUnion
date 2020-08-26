#pragma once

#include "BaseDecoder.h"
#include "Domain.h"

class ScDecoder : public BaseDecoder {

private:
public:
	ScDecoder(PolarCode * code);

	std::vector<int> Decode(std::vector<double> llr) override;
	domain GetDomain() override;

	~ScDecoder() {};
};