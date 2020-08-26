#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0
 

ScDecoder::ScDecoder(PolarCode * codePtr) :  BaseDecoder(codePtr) {
}

domain ScDecoder::GetDomain() {
	return LLR;
}

double f(double x, double y) {
	double sign = 1.0;

	if (x < 0) {
		sign *= -1;
		x *= -1;
	}
	if (y < 0) {
		sign *= -1;
		y *= -1;
	}
		
	return ((x < y) ? x  : y ) * sign;
}

double g(double x, double y, int b) {
	return y + (1 - 2 * b) * x;
}


std::vector<double>  GetLeftMessage(std::vector<double>::iterator leftIt, std::vector<double>::iterator rightIt, size_t n) {

	std::vector<double> out(n, 0);
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++)
	{
		out[i] = f(*leftIt, *rightIt);
	}

	return out;
}

std::vector<double> GetRightMessage(std::vector<double>::iterator leftIt, std::vector<double>::iterator rightIt, std::vector<int> uHat, size_t n) {
	std::vector<double> out(n, 0);

	for (size_t i = 0; i < n; i++, leftIt++, rightIt++)
	{
		out[i] = g(*leftIt, *rightIt, uHat[i]);
	}

	return out;
}

// return u, cause they are needed at each step
std::vector<int> DecodeRecursive(
	std::vector<double>::iterator inLlrIt,
	std::vector<double>::iterator outLlrIt,
	std::vector<int>::iterator maskIt,
	std::vector<int>::iterator xIt,
	size_t n
	) 
{
	std::vector<int> uHat(n, 0);
	if (n == 1) {
		*outLlrIt = *inLlrIt;
		if (*maskIt)
			*xIt = *inLlrIt < 0;

		else
			*xIt = FROZEN_VALUE;
		
		uHat[0] = *xIt;
	}
	else {
		size_t half_n = n / 2;
		std::vector<double> llrLeft = GetLeftMessage(inLlrIt, inLlrIt + half_n, half_n);
		std::vector<int> uHatLeft = DecodeRecursive(llrLeft.begin(), outLlrIt, maskIt, xIt, half_n);

		std::vector<double> llrRight = GetRightMessage(inLlrIt, inLlrIt + half_n, uHatLeft, half_n);
		std::vector<int> uHatRight = DecodeRecursive(llrRight.begin(), outLlrIt + half_n, maskIt + half_n, xIt + half_n, half_n);

		for (size_t i = 0; i < half_n; i++)
		{
			uHat[i] = (uHatLeft[i] + uHatRight[i]) % 2;
			uHat[i + half_n] = uHatRight[i];
		}
		
	}
	return uHat;
}

std::vector<int> ScDecoder::Decode(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	std::vector<double> outLlr(n, 0);
	std::vector<int> x(n, 0);
	std::vector<int> mask = _codePtr->BitsMask();

	std::vector<int> u = DecodeRecursive(inLlr.begin(), outLlr.begin(), mask.begin(), x.begin(), n);

	size_t k = _codePtr->k();
	std::vector<int> result(k, 0);
	size_t j = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (mask[i] == 0)
			continue;

		result[j] = x[i];
		j++;
	}

	return result;
}