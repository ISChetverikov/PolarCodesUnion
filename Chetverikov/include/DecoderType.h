#pragma once

#include <string>
#include <unordered_map>

enum decoderType { SC, SCFano, UnknownDecoder, SCFlip };

decoderType decoderTypeFromString(std::string str);
std::string decoderTypeToString(decoderType type);