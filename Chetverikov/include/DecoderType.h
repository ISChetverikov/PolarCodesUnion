#pragma once

#include <string>
#include <unordered_map>

enum decoderType { SC, SCRecursive, SCFano, UnknownDecoder, SCFlip, SCFlipProg, SCList, SCListFlipStat, SCListFlipOracleStat, SCFlipFano };

decoderType decoderTypeFromString(std::string str);
std::string decoderTypeToString(decoderType type);