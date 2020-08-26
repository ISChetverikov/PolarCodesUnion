#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScFlipDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0


ScFlipDecoder::ScFlipDecoder(PolarCode * codePtr, int T) : BaseDecoder(codePtr) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	_treeHeight = m + 1;

	for (size_t i = 0; i < _treeHeight; i++)
	{
		std::vector<double> b(n, -10000.0);
		std::vector<int> u(n, -1);

		_beliefTree.push_back(b);
		_uhatTree.push_back(u);
	}

	_mask = _codePtr->BitsMask();
	_x = std::vector<int>(n, -1);
	_T = T;
}

domain ScFlipDecoder::GetDomain() {
	return LLR;
}

double ScFlipDecoder::f(double llr1, double llr2) {
	double sign = 1.0;

	if (llr1 < 0) {
		sign *= -1;
		llr1 *= -1;
	}
	if (llr2 < 0) {
		sign *= -1;
		llr2 *= -1;
	}

	return ((llr1 < llr2) ? llr1 : llr2) * sign;
}

double ScFlipDecoder::g(double llr1, double llr2, int u1) {
	return llr2 + (1 - 2 * u1) * llr1;
}

int ScFlipDecoder::L(double llr) {
	return (llr >= 0) ? 0 : 1;
}

void ScFlipDecoder::FillLeftMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++)
	{
		*outIt = f(*leftIt, *rightIt);
	}
}

void ScFlipDecoder::FillRightMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<int>::iterator uhatIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++, uhatIt++)
	{
		*outIt = g(*leftIt, *rightIt, *uhatIt);
	}
}


size_t ScFlipDecoder::log2(size_t n) {
	size_t m = 0;
	while (n > 0)
	{
		n = n >> 1;
		m++;
	}

	return m;
}

void ScFlipDecoder::PassDown(size_t iter) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();

	size_t iterXor;
	size_t level;
	if (iter) {
		iterXor = iter ^ (iter - 1);
		level = m - log2(iterXor);
	}
	else {
		level = 0;
	}

	std::vector<int> binaryIter(m - level, 0);
	size_t iterCopy = iter;
	for (int i = (int)binaryIter.size() - 1; i >= 0; i--)
	{
		binaryIter[i] = iterCopy % 2;
		iterCopy = iterCopy >> 1;
	}

	size_t length = (size_t)1 << (m - level - 1);
	for (size_t i = level; i < m; i++)
	{
		size_t ones = ~0u;
		size_t offset = iter & (ones << (m - i));

		if (!binaryIter[i - level]) {
			FillLeftMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_beliefTree[i + 1].begin() + offset,
				length);
		}
		else {

			FillRightMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_uhatTree[i + 1].begin() + offset,
				_beliefTree[i + 1].begin() + offset + length,
				length);
		}

		length = length / 2;
	}
}

void ScFlipDecoder::PassUp(size_t iter) {
	size_t m = _codePtr->m();
	size_t iterCopy = iter;

	size_t bit = iterCopy % 2;
	size_t length = 1;
	size_t level = m;
	size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + i] = _uhatTree[level][offset + i] ^ _uhatTree[level][offset + length + i];
		}
		for (size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + length + i] = _uhatTree[level][offset + length + i];
		}


		iterCopy = iterCopy >> 1;
		bit = iterCopy % 2;
		length *= 2;
		level -= 1;
	}
}

std::vector<int> ScFlipDecoder::Decode(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();
	size_t k = _codePtr->k();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inLlr[i];
	}

	for (size_t i = 0; i < n; i++)
	{
		PassDown(i);
		if (_mask[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}

	std::vector<int> result(k, 0);
	size_t l = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_mask[i] == 0)
			continue;

		result[l] = _x[i];
		l++;
	}

	return result;
}