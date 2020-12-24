
#include <random>
#include "../include/BpskAwgnChannel.h"

double snrToSigma(double snr) {
	return sqrt(pow(10, -snr / 10) / 2);
}

BpskAwgnChannel::BpskAwgnChannel() : BaseChannel() {
	_snr = 1.0; // trash value
	_sigma = snrToSigma(_snr);

	std::normal_distribution<double> normal_dist(0, _sigma);
	_normal_dist = normal_dist;
}



void BpskAwgnChannel::SetSnr(double snr) {
	BaseChannel::SetSnr(snr);
	_sigma = snrToSigma(snr);
	std::normal_distribution<double> normal_dist(0, _sigma);
	_normal_dist = normal_dist;
}

double LlrToP1(double llr) {
	if (llr > 300.0)
		return 0.0;

	if (llr < -300.0)
		return 1.0;

	return 1.0 / (1 + exp(llr));
}

double InputToLlr(double input, double sigma) {
	return 2 * input / (sigma * sigma);
}

int ModulateBpsk(int input) {
	return 1 - 2 * input;
}

std::vector<double> BpskAwgnChannel::Pass(std::vector<int> codeword) {
	size_t n = codeword.size();
	std::vector<double> output(n, 0);
	for (size_t i = 0; i < n; i++) {
		output[i] = ModulateBpsk(codeword[i]) + _normal_dist(_randomDevice);
		output[i] = InputToLlr(output[i], _sigma);
	}

#ifdef DOMAIN_P1

	for (size_t i = 0; i < n; i++)
	{
		output[i] = LlrToP1(output[i]);
	}
	
#endif // DOMAIN_P1

	return output;

}