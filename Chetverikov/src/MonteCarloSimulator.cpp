#include <random>
#include <chrono>
#include <exception>
#include "../include/SimulationParameters.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/ScListDecoder.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTestsCount,
	int maxRejectionsCount,
	PolarCode * codePtr,
	Encoder * encoderPtr,
	BaseDecoder * decoderPtr,
	bool isSigmaDependOnR) : BaseSimulator(codePtr, encoderPtr, decoderPtr, isSigmaDependOnR)
{
	_maxTestsCount = maxTestsCount;
	_maxRejectionsCount = maxRejectionsCount;
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

SimulationIterationResults MonteCarloSimulator::Run(double snr)
{	
	SimulationIterationResults result;

	std::random_device randomDevice;
	size_t n = _codePtr->N();
	size_t k = _codePtr->k();
	std::vector<int> word(k, 0);
	std::vector<int> codeword(n, 0);
	std::vector<double> beliefs(n, 0);
	std::vector<int> decoded(n, 0);

	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr, (double)k / n);
	int tests = 0;
	int wrong_dec = 0;
	std::normal_distribution<double> normal_dist(0, sigma);
	std::uniform_int_distribution<> uniform_discrete_dist(0, 1);

	_decoderPtr->SetSigma(sigma);

	while ((tests < _maxTestsCount || _maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount)) {
		tests++;

		std::generate(word.begin(), word.end(), [&]() { return uniform_discrete_dist(randomDevice); });

		codeword = _encoderPtr->Encode(word);
		// Give answer to a decoder for debugging or statistic retreving
		_decoderPtr->SetCodeword(_encoderPtr->PolarTransform(codeword));

		for (size_t i = 0; i < n; i++) {
			beliefs[i] = ModulateBpsk(codeword[i]) + normal_dist(randomDevice);
			beliefs[i] = InputToLlr(beliefs[i], sigma);
		}

		auto domain = _decoderPtr->GetDomain();
		if (domain == P1) {
			for (size_t i = 0; i < n; i++)
			{
				beliefs[i] = LlrToP1(beliefs[i]);
			}
		}
		
		decoded = _decoderPtr->Decode(beliefs);
		if (tests % 1000 == 0)
			std::cout << tests << std::endl;
		// LIST4 parallel decoder
		/*auto d = new ScListDecoder(_codePtr, 4);
		auto de = d->Decode(beliefs);

		if (false && de != decoded && de == word) {
			std::cout << "Belief:\n{";
			for (size_t i = 0; i < beliefs.size(); i++)
			{
				std::cout << beliefs[i] << ", ";
			}
			std::cout << "}" << std::endl;

			std::cout << "Codeword:\n{";
			for (size_t i = 0; i < codeword.size(); i++)
			{
				std::cout << codeword[i] << ", ";
			}
			std::cout << "}" << std::endl;

			std::vector<int> extWord = _encoderPtr->PolarTransform(codeword);

			std::cout << "ExtWord:\n{";
			for (size_t i = 0; i < extWord.size(); i++)
			{
				std::cout << extWord[i] << ", ";
			}
			std::cout << "}" << std::endl;

			std::cout << "Word:\n{";
			for (size_t i = 0; i < word.size(); i++)
			{
				std::cout << word[i] << ", ";
			}
			std::cout << "}\n";

			throw "Exit debug!";
		}*/

		if (decoded != word)
			wrong_dec += 1;
	}

	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, k, n);
	result.sigma = sigma;

	result.fer = (double)wrong_dec / tests;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = tests;
	result.rejectionsCount = wrong_dec;
	
	return result;
}
