#include <math.h>

#include "../include/PolarCode.h"
#include "../include/Exceptions.h"
#include "../include/Encoder.h"
#include "../include/CRC.h"

Encoder::Encoder(PolarCode * codePtr) {
	_codePtr = codePtr;
}

void PolarTransform(std::vector<int>::iterator begin, size_t n) {
	
	if (n == 1)
		return;

	size_t half_n = n / 2;
	for (auto it = begin; it < begin + half_n; it++)
		*it = (*it + *(it + half_n)) % 2;
	
	PolarTransform(begin, half_n);
	PolarTransform(begin + half_n, half_n);

	return;
}

std::vector<int> Encoder::Encode(std::vector<int> word) {
	size_t k = word.size();

	if (k != _codePtr->k())
		throw IncorrectVectorSizeException("Length of the word does not matched with code parameter k");

	size_t n = _codePtr->N();
	std::vector<int> codeword(n, 0);
	std::vector<int> mask = _codePtr->BitsMask();

	size_t j = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (mask[i] == 0)
			continue;

		codeword[i] = word[j];
		j++;
	}

	PolarTransform(codeword.begin(), n);
	return codeword;
}