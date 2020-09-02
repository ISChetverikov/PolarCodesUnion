#pragma once
#include <vector>
class PolarCode {

protected:
	
	size_t _m = 0;
	size_t _N = 0;
	size_t _k = 0;
	// with length of CRC hash
	size_t _k_extended = 0;

	std::vector<int> _bitsMask;
	std::vector<int> _crcMask;

	std::vector<int> _unfrozenBits;
	std::vector<int> _crcUnfrozenBits;

	std::vector<int> _crcPoly;
	size_t _crcDeg;

public:
	PolarCode();
	PolarCode(int m, int k, std::vector<int> reliabilitySequence, std::vector<int> crcPoly);

	size_t m();
	size_t N();
	size_t k();
	
	std::vector<int> BitsMask();
	std::vector<int> CrcMask();

	std::vector<int> UnfrozenBits();
	std::vector<int> CrcUnfrozenBits();

	std::vector<int> CrcPoly();
	size_t CrcDeg();

	bool IsCrcUsed();

	~PolarCode() {};
};