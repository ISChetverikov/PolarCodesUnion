#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScFlipDecoder.h"
#include "../include/PolarCode.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define FROZEN_VALUE 0

ScFlipDecoder::ScFlipDecoder(PolarCode * codePtr, int T) : ScCrcAidedDecoder(codePtr) {
	_T = T;
}

std::vector<int> ScFlipDecoder::GetSmallestBeliefsIndices(std::vector<double> beliefs, int count) {
	std::vector<int> indices(count, 0);
	for (size_t i = 0; i < beliefs.size(); i++)
	{
#ifdef DOMAIN_LLR
		beliefs[i] = fabs(beliefs[i]);
#elif DOMAIN_P1
		beliefs[i] = fabs(beliefs[i] - 0.5);
#endif // DOMAIN_LLR

	}

	for (size_t i = 0; i < count; i++)
	{
		auto minIt = std::min_element(beliefs.begin(), beliefs.end());
		auto minInd = (int)std::distance( beliefs.begin(), minIt);
		indices[i] = minInd;
		beliefs.erase(minIt);
	}
	
	return indices;
}

std::vector<int> ScFlipDecoder::Decode(std::vector<double> inBeliefs) {

	DecodeInternal(inBeliefs);

	size_t m = _codePtr->m();
	if (_T > 0 && !IsCrcPassed(_x)) {
		std::vector<int> suspectedBits = GetSmallestBeliefsIndices(_beliefTree[m], _T);
		for (size_t i = 0; i < _T; i++)
		{
			int bitPosition = suspectedBits[i];
			_x[bitPosition] = !_x[bitPosition];

			DecodeFrom(bitPosition);

			if (IsCrcPassed(_x))
				break;

			_x[bitPosition] = !_x[bitPosition];
			DecodeFrom(bitPosition);
		}
	}

	return TakeResult();
}