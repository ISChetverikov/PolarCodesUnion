
#include <random>
#include "../include/BscChannel.h"

BscChannel::BscChannel(double coderate) : BaseChannel() {
	std::bernoulli_distribution bernoulli_dist(0.5);
	_bernoulli_dist = bernoulli_dist;
	
	_coderate = coderate;
	_fixedLlr = 1.0;
}

double snrToEbN0(double snr, double coderate) {
	return snr - 10 * log10(coderate);
}

double ebnoToPErr(double ebno) {
	return 1.0 / 2 * erfc(sqrt(ebno));
}


void BscChannel::SetSnr(double snr) {
	BaseChannel::SetSnr(snr);

	double p = ebnoToPErr(snrToEbN0(_snr, _coderate));
	_fixedLlr = log((1 - p) / p);
	std::bernoulli_distribution b(p);
	_bernoulli_dist = b;
}

double LlrToP12(double llr) {
	if (llr > 300.0)
		return 0.0;

	if (llr < -300.0)
		return 1.0;

	return 1.0 / (1 + exp(llr));
}


std::vector<double> BscChannel::Pass(std::vector<int> codeword) {
	size_t n = codeword.size();
	std::vector<double> output(n, 0);

	for (size_t i = 0; i < n; i++) {

		if (_bernoulli_dist(_randomDevice)) {
			output[i] = !codeword[i];
		}
		else
		{
			output[i] = codeword[i];
		}

		if (output[i] == 1) {
			output[i] = -_fixedLlr;
		}
		else {
			output[i] = _fixedLlr;
		}
	}

#ifdef DOMAIN_P1

	for (size_t i = 0; i < n; i++)
	{
		output[i] = LlrToP12(output[i]);
	}

#endif // DOMAIN_P1

	return output;

}