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
	_crcDeg = _crcPoly.size();

	if (_k + _crcDeg > _N)
		throw CrcPolyException("Dimensions of a message with crc do not fit with dimension of a codeword ");

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
	if (IsCrcUsed()) {
		_crcMask = std::vector<int>(_N, 0);
		_crcUnfrozenBits = std::vector<int>(_crcDeg, 0);
		for (size_t i = sequenceLength - _k_extended, j = 0; i < sequenceLength - k; i++, j++)
		{
			_crcMask[reliabilitySequence[i]] = 1;
			_crcUnfrozenBits[j] = reliabilitySequence[i];
		}
	}
	
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

std::vector<int> PolarCode::CrcMask() {
	return _crcMask;
}

std::vector<int> PolarCode::CrcUnfrozenBits() {
	return _crcUnfrozenBits;
}

bool PolarCode::IsCrcUsed() {
	return !_crcPoly.empty();
}

std::vector<int> PolarCode::CrcPoly() {
	return _crcPoly;
}

size_t PolarCode::CrcDeg() {
	return _crcDeg;
}