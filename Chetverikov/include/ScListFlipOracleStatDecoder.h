#pragma once

#include <map>
#include "ScListFlipStatDecoder.h"

class ScListFlipOracleStatDecoder : public ScListFlipStatDecoder {
private:
	std::map<int, int> _singleOracleFlipsStat;
	std::map<std::tuple<int, int>, int> _doubleOracleFlipsStat;
	std::map<std::tuple<int, int, int>, int> _tripleOracleFlipsStat;
	
protected:
	int GetFirstErrorPosition(std::vector<int> codeword1, std::vector<int> codeword2);

public:
	ScListFlipOracleStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipOracleStatDecoder() {};
};