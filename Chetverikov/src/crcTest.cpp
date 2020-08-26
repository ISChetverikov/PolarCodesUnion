#include "../include/CRC.h"
#include <vector>
#include <iostream>

int main() {
	std::vector<int> poly = { 1, 1, 1, 0, 0, 0, 0, 0, 1};
	//std::vector<int> input = { 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<int> input = { 1, 1, 1, 0};

	auto c = CRC(poly);
	std::vector<int> hash = c.Calculate(input);
	for (size_t i = 0; i < hash.size(); i++)
	{
		std::cout << hash[i];
	}
	std::cout << std::endl;

	return 0;
}