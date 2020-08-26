#pragma once 

#include <vector>

class CRC {
private:
	std::vector<int> _poly;
	std::vector<int> _init;

	size_t _deg = 0;
	int _paddingSymbol = 0;

public:

	CRC(std::vector<int> poly);
	std::vector<int> Calculate(std::vector<int> bits);

	std::vector<int> Add(std::vector<int> bits);
	bool Check(std::vector<int> bits);
};