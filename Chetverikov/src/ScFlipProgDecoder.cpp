#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/GaussianApproximation.h"
#include "../include/ScFlipProgDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0
#define ROOT_POSITION -1

// mask: 1 - black, 0 - white
std::vector<int> ScFlipProgDecoder::GetCriticalSet(std::vector<int> mask, int position) {
	std::vector<std::vector<int>> tree;
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	std::vector<int> result;

	for (size_t i = 0; i < m + 1; i++)
	{
		std::vector<int> temp((int)1 << i, 0);
		tree.push_back(temp);
	}
	for (size_t i = position; i < n; i++)
	{
		tree[m][i] = mask[i];
	}

	for (int i = (int)m - 1; i >= 0; i--)
	{
		size_t levelSize = tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			int left = j << 1;
			int right = left + 1;

			tree[i][j] = tree[i + 1][left] && tree[i + 1][right];
		}
	}

	for (size_t i = 0; i < m + 1; i++)
	{
		int levelSize = (int)tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			if (!tree[i][j])
				continue;

			int criticalBit = j << (m - i);
			result.push_back(criticalBit);

			int left = j;
			int right = j;
			tree[i][j] = 0;
			for (size_t k = i + 1; k <= m; k++)
			{
				left = left << 1;
				right = (right << 1) + 1;

				for (size_t l = left; l <= right; l++)
					tree[k][l] = 0;
			}
		}
	}

	return result;
}

CriticalSetNode * ScFlipProgDecoder::GetCriticalSetTree(std::vector<int> mask, int levelMax) {
	std::queue<CriticalSetNode *> q;

	auto root = new CriticalSetNode();
	root->Bit = ROOT_POSITION;
	q.push(root);
	while (!q.empty()) {
		CriticalSetNode * cur = q.front();
		q.pop();

		auto criticalSet = GetCriticalSet(_maskWithCrc, cur->Bit + 1);
		for (size_t i = 0; i < criticalSet.size(); i++)
		{
			auto childNode = new CriticalSetNode();
			childNode->Bit = criticalSet[i];
			
			if (cur->Bit != ROOT_POSITION) {
				childNode->Path = cur->Path;
				childNode->Path.push_back(cur->Bit);
			}
				
			cur->Children.push_back(childNode);
			
			if (childNode->Path.size() < levelMax - 1) // indices of level from 1
				q.push(childNode);
		}
	}

	return root;
}

ScFlipProgDecoder::ScFlipProgDecoder(PolarCode * codePtr, int level, double gammaLeft, double gammaRight, std::vector<double> omegaArr) : BaseDecoder(codePtr) {
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

	_maskWithCrc = std::vector<int>(n, 0);
	auto maskInf = _codePtr->BitsMask();
	
	if (_codePtr->IsCrcUsed()) {
		auto maskCrc = _codePtr->CrcMask();
		for (size_t i = 0; i < n; i++)
			_maskWithCrc[i] = maskInf[i] || maskCrc[i];
	}	
	else
		_maskWithCrc = maskInf;
	
	_gammaLeft = gammaLeft;
	_gammaRight = gammaRight;
	_levelMax = level;
	_omegaArr = omegaArr;

	_criticalSetTree = GetCriticalSetTree(_maskWithCrc, _levelMax);
	_x = std::vector<int>(n, -1);
	_crcPtr = new CRC(_codePtr->CrcPoly());
	_subchannelsMeansGa = std::vector<double>(n, 0);

	// optimization allocation
	_binaryIter = std::vector<int>(m, 0);
}

domain ScFlipProgDecoder::GetDomain() {
	return LLR;

}

void ScFlipProgDecoder::SetSigma(double sigma) {
	BaseDecoder::SetSigma(sigma);

	GaussianApproximation ga(sigma);
	size_t n = _codePtr->N();
	for (size_t i = 0; i < n; i++)
	{
		_subchannelsMeansGa[i] = ga.GetMu(i + 1, n);
	}

	return;
}

//double ScFlipProgDecoder::f(double llr1, double llr2) {
//	double prod = tanh(llr1 / 2) * tanh(llr2 / 2);
//	double limit = 0.9999999999999999;
//
//	if (prod > limit)
//		prod = limit;
//	if (prod < -limit)
//		prod = -limit;
//
//	return 2 * atanh(prod);
//}


double ScFlipProgDecoder::f(double llr1, double llr2) {
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

double ScFlipProgDecoder::g(double llr1, double llr2, int u1) {
	return llr2 + (1 - 2 * u1) * llr1;
}

int ScFlipProgDecoder::L(double llr) {
	return (llr >= 0) ? 0 : 1;
}

void ScFlipProgDecoder::FillLeftMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++)
	{
		*outIt = f(*leftIt, *rightIt);
	}
}

void ScFlipProgDecoder::FillRightMessageInTree(std::vector<double>::iterator leftIt,
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


size_t ScFlipProgDecoder::log2(size_t n) {
	size_t m = 0;
	while (n > 0)
	{
		n = n >> 1;
		m++;
	}

	return m;
}

void ScFlipProgDecoder::PassDown(size_t iter) {
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

	//std::vector<int> binaryIter(m - level, 0);
	int size = m - level;
	size_t iterCopy = iter;
	for (int i = size - 1; i >= 0; i--)
	{
		_binaryIter[i] = iterCopy % 2;
		iterCopy = iterCopy >> 1;
	}

	size_t length = (size_t)1 << (m - level - 1);
	for (size_t i = level; i < m; i++)
	{
		size_t ones = ~0u;
		size_t offset = iter & (ones << (m - i));

		if (!_binaryIter[i - level]) {
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

void ScFlipProgDecoder::PassUp(size_t iter) {
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

void ScFlipProgDecoder::DecodeFromTo(int startPosition, int endPosition) {
	size_t n = _codePtr->N();
	size_t m = _codePtr->m();

	_uhatTree[m][startPosition] = _x[startPosition];
	PassUp(startPosition);
	for (size_t i = startPosition + 1; i < endPosition; i++)
	{
		PassDown(i);
		if (_maskWithCrc[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}

	if (endPosition != n) {
		PassDown(endPosition);
		if (_maskWithCrc[endPosition]) {
			_x[endPosition] = L(_beliefTree[m][endPosition]);
		}
		else {
			_x[endPosition] = FROZEN_VALUE;
		}
	}
		
}

bool ScFlipProgDecoder::IsCrcPassed(std::vector<int> codeword) {
	size_t n = codeword.size();
	size_t k = _codePtr->k();
	size_t deg = _codePtr->CrcDeg();

	auto wordBits = _codePtr->UnfrozenBits();
	std::vector<int> word(k, 0);
	for (size_t i = 0; i < k; i++)
	{
		word[i] = codeword[wordBits[i]];
	}

	auto crcBits = _codePtr->CrcUnfrozenBits();
	std::vector<int> crc(deg, 0);
	for (size_t i = 0; i < deg; i++)
	{
		crc[i] = codeword[crcBits[i]];
	}

	auto crcReal = _crcPtr->Calculate(word);
	return crc == crcReal;
}


// with using means after gaussian approxiamtion procedure
std::vector<CriticalSetNode *> ScFlipProgDecoder::SortCriticalNodes(std::vector<CriticalSetNode *> criticalNodes, std::vector<double> llrs) {
	size_t length = criticalNodes.size();
	std::vector<CriticalSetNode *> result(length, 0);
	std::vector<double> sortingMetric(length, 0);
	
	for (size_t i = 0; i < length; i++)
	{
		int criticalInd = criticalNodes[i]->Bit;
		sortingMetric[i] = fabs(llrs[criticalInd]) / _subchannelsMeansGa[criticalInd];
	}

	for (size_t i = 0; i < length; i++)
	{
		auto minIt = std::min_element(sortingMetric.begin(), sortingMetric.end());
		int minInd = (int)std::distance(sortingMetric.begin(), minIt);
		result[i] = criticalNodes[minInd];

		criticalNodes.erase(criticalNodes.begin() + minInd);
		sortingMetric.erase(sortingMetric.begin() + minInd);
	}

	return result;
}

bool ScFlipProgDecoder::NoChild(CriticalSetNode * node, std::vector<double> inLlr) {
	int position = node->Bit;
	size_t n = _codePtr->N();
	size_t criticalSetSize = node->Children.size();
	
	if (position >= n)
		return true;

	if (criticalSetSize <= 0)
		return true;

	size_t level = node->Path.size() + 1;

	if (_omegaArr[level] >= 1.0)
		return false;

	std::vector<int> maskWithoutCriticalSet = _maskWithCrc;
	for (size_t i = 0; i < criticalSetSize; i++)
	{
		maskWithoutCriticalSet[node->Children[i]->Bit] = 0;
	}
	maskWithoutCriticalSet[node->Bit] = 0;

	int n1 = 0;
	int n2 = 0;
	for (size_t i = position; i < n; i++)
	{
		if (maskWithoutCriticalSet[i]) {
			n1++;

			double mu = _subchannelsMeansGa[i];
			double sigma = sqrt(2 * mu);
			if (inLlr[i] < mu - _gammaLeft * sigma)
				n2++;
		}
	}

	return ((double)n2 ) / n1 > _omegaArr[level];
}
bool ScFlipProgDecoder::NotSelect(CriticalSetNode * node, std::vector<double> inLlr) {
	int position = node->Bit;
	size_t level_dec = node->Path.size();

	if (_omegaArr[level_dec] >= 1.0)
		return false;

	//std::cout << level_dec << std::endl;

	double mu = _subchannelsMeansGa[position];
	double sigma = sqrt(2 * mu);
	
	return inLlr[position] > mu + _gammaRight * sigma;
}

std::vector<int> ScFlipProgDecoder::Decode(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();
	
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inLlr[i];
	}

	for (size_t i = 0; i < n; i++)
	{
		PassDown(i);
		if (_maskWithCrc[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}

	if (!IsCrcPassed(_x)) {
		std::queue<CriticalSetNode *> q;
		q.push(_criticalSetTree);

		bool exitFlag = false;
		while (!q.empty() && !exitFlag) {
			auto currentNode = q.front();
			q.pop();

			auto criticalNodes = currentNode->Children;
			std::vector<CriticalSetNode *> suspectedNodes = SortCriticalNodes(criticalNodes, _beliefTree[m]);

			for (size_t i = 0; i < suspectedNodes.size(); i++)
			{
				int bitPosition = suspectedNodes[i]->Bit;
				if (NotSelect(suspectedNodes[i], inLlr))
					continue;

				size_t pathSize = suspectedNodes[i]->Path.size();
				for (size_t j = 0; j < pathSize; j++)
				{
					int prevBitPosition = suspectedNodes[i]->Path[j];
					int nextPrevBitPosition = (j != pathSize - 1) ? suspectedNodes[i]->Path[j + 1] : bitPosition;
					_x[prevBitPosition] = !_x[prevBitPosition];
					DecodeFromTo(prevBitPosition, nextPrevBitPosition);
				}
				
				_x[bitPosition] = !_x[bitPosition];
				DecodeFromTo(bitPosition, (int)n);
				
				if (IsCrcPassed(_x)) {
					exitFlag = true;
					break;
				}
				
				int minPosition = suspectedNodes[i]->Path.empty() ? bitPosition : suspectedNodes[i]->Path[0];
				_x[minPosition] = !_x[minPosition];
				DecodeFromTo(minPosition, (int)n);
				
				if (!NoChild(suspectedNodes[i], inLlr))
					q.push(suspectedNodes[i]);
			}
		}

	}

	std::vector<int> result(_codePtr->k(), 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _x[codewordBits[i]];
	}

	return result;
}