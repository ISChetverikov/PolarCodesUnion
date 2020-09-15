#include "../include/CRC.h"
#include <vector>
#include <iostream>
#include <cmath>

double f2(double llr1, double llr2) {
	double prod = tanh(llr1 / 2) * tanh(llr2 / 2);
	double limit = 0.9999999999999999;

	if (prod > limit)
		prod = limit;
	if (prod < -limit)
		prod = -limit;


	return 2 * atanh(prod);
}

int main1() {
	std::vector<int> poly = { 1, 1, 1};
	//std::vector<int> input = { 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<int> input = { 1, 0, 0, 0 };
	auto c = CRC(poly);
	std::vector<int> hash = c.Calculate(input);

	int length = input.size();
	for (int i = 0; i < length; i++)
	{
		
		for (int j = 0; j < length; j++)
		{
			
			std::vector<int> newInput = input;
			newInput[i] = !input[i];
			newInput[j] = !input[j];

			std::vector<int> newHash = c.Calculate(newInput);

			std::cout << (newHash == hash);
		}
		std::cout << std::endl;
	}


	return 0;
}