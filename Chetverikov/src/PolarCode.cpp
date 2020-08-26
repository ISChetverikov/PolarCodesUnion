#include "../include/PolarCode.h"
#include "../include/Exceptions.h"
#include <algorithm>
PolarCode::PolarCode() {

}
PolarCode::PolarCode(int m, int k, std::vector<int> reliabilitySequence, std::vector<int> crcPoly) {
	_m = m;
	_k = k;
	_N = 2 << (m - 1);
	size_t sequenceLength = reliabilitySequence.size();
	if (_N != sequenceLength)
		throw IncorrectSequenceSizeException("Sequence length does not match with code length");

	_crcPoly = crcPoly;
	_crcDeg = _crcPoly.size() - 1;

	if (!_crcPoly.empty())
		_k_extended = _k + _crcDeg;

	// the mask and the set without crc
	_bitsMask = std::vector<int>(_N, 0);
	_unfrozenBits = std::vector<int>(k, 0);
	for (size_t i = sequenceLength - k, j = 0; i < sequenceLength; i++, j++)
	{
		_bitsMask[reliabilitySequence[i]] = 1;
		_unfrozenBits[j] = reliabilitySequence[i];
	}
	sort(_unfrozenBits.begin(), _unfrozenBits.end());

	// the mask and the set with crc
	_bitsMaskExtended = std::vector<int>(_N, 0);
	_unfrozenBitsExtended = std::vector<int>(_k_extended, 0);
	for (size_t i = sequenceLength - _k_extended, j = 0; i < sequenceLength; i++, j++)
	{
		_bitsMaskExtended[reliabilitySequence[i]] = 1;
		_unfrozenBitsExtended[j] = reliabilitySequence[i];
	}
	sort(_unfrozenBitsExtended.begin(), _unfrozenBitsExtended.end());

}

size_t PolarCode::m() {
	return _m;
}
size_t PolarCode::N() {
	return _N;
}
size_t PolarCode::k() {
	return _k;
}
std::vector<int> PolarCode::BitsMask() {
	return _bitsMask;
}

std::vector<int> PolarCode::UnfrozenBits() {
	return _unfrozenBits;
}

std::vector<int> PolarCode::BitsMaskExtended() {
	return _bitsMaskExtended;
}

std::vector<int> PolarCode::UnfrozenBitsExtended() {
	return _unfrozenBitsExtended;
}

bool PolarCode::IsCrcUsed() {
	return _crcPoly.empty();
}

std::vector<int> PolarCode::CrcPoly() {
	return _crcPoly;
}