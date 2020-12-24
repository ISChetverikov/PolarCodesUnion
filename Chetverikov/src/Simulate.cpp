#include <exception>
#include <chrono>
#include <string>
#include <fstream>

#include "../include/Simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/BaseDecoder.h"

#include "../include/MonteCarloSimulator.h"
#include "../include/PolarCode.h"
#include "../include/Encoder.h"
#include "../include/ScRecursiveDecoder.h"
#include "../include/ScDecoder.h"
#include "../include/ScFanoDecoder.h"
#include "../include/ScFlipDecoder.h"
#include "../include/ScListFanoDecoder.h"
#include "../include/ScListFlipStatDecoder.h"
#include "../include/ScListFlipOracleStatDecoder.h"
#include "../include/ScListDecoder.h"
#include "../include/ScFlipProgDecoder.h"
#include "../include/BaseSimulator.h"
#include "../include/ConfigReading.h"
#include "../include/Exceptions.h"
#include "../include/BpskAwgnChannel.h"
#include "../include/BscChannel.h"

int ExtractInt(std::unordered_map<std::string, std::string> map, std::string key, std::string section) {
	if (map.count(key) <= 0)
		throw MissedParamException("Missed parameters \"" + key + "\" in section " + section);

	int result;
	try {
		result = std::stoi(map[key]);
	}
	catch (std::exception& e) {
		throw ParseParamException("Parse error of parameter \"" + key + "\": " + e.what());
	}
	
	return result;
}

double ExtractDouble(std::unordered_map<std::string, std::string> map, std::string key, std::string section) {
	if (map.count(key) <= 0)
		throw MissedParamException("Missed parameters \"" + key + "\" in section " + section);

	double result;
	try {
		result = std::stod(map[key]);
	}
	catch (std::exception& e) {
		throw ParseParamException("Parse error of parameter \"" + key + "\": " + e.what());
	}

	return result;
}

std::string ExtractString(std::unordered_map<std::string, std::string> map, std::string key, std::string section, bool isRequired) {
	if (map.count(key) <= 0) {
		if (isRequired)
			throw MissedParamException("Missed parameters \"" + key + "\" in section " + section);

		else
			return "";
	}
		

	std::string result;
	try {
		result = map[key];
	}
	catch (std::exception& e) {
		throw ParseParamException("Parse error of parameter \"" + key + "\": " + e.what());
	}

	return result;
}

std::vector<int> ReadSequenceFromFile(std::string path) {
	std::vector<int> seq;
	std::string line;
	std::ifstream myFile(path);

	std::getline(myFile, line);

	int val;
	std::stringstream ss(line);

	while (ss >> val)
		seq.push_back(val);

	return seq;
}

std::vector<int> PolyStrToVector(std::string crcPolyStr) {
	if (crcPolyStr.empty())
		return std::vector<int>();

	size_t deg = crcPolyStr.size();

	std::vector<int> result(deg, 0);

	for (size_t i = 0; i < deg; i++)
	{
		if (crcPolyStr[i] == '1')
			result[deg - i - 1] = 1;

		else if (crcPolyStr[i] != '0')
			throw CrcPolyException("Crc polynom has incorrect symbol: " + crcPolyStr);
	}

	return result;
}

std::vector<double> OmegaArrStrToVector(std::string str)
{
	if (str.empty())
		return std::vector<double>();

	std::vector<double> result;
	std::string token;
	std::istringstream tokenStream(str);
	while (std::getline(tokenStream, token, ','))
	{
		result.push_back(std::stod(token));
	}
	return result;
}

PolarCode * BuildCode(std::unordered_map<std::string, std::string> codeParams) {
	PolarCode * codePtr;

	int m = ExtractInt(codeParams, "m", "PolarCode");
	int k = ExtractInt(codeParams, "k", "PolarCode");
	std::string crcPolyString = ExtractString(codeParams, "CRC", "PolarCode", false);
	std::vector<int> crcPoly = PolyStrToVector(crcPolyString);
	std::string sequenceFilePath = ExtractString(codeParams, "sequenceFile", "PolarCode", true);
	std::vector<int> reliabilitySequence = ReadSequenceFromFile(sequenceFilePath);

	codePtr = new PolarCode(m, k, reliabilitySequence, crcPoly);
	return codePtr;
}

BaseDecoder * BuildDecoder(
                            decoderType decoderType,
                            std::unordered_map<std::string, std::string> decoderParams,
							PolarCode * codePtr) {
    BaseDecoder * decoderPtr = NULL;
	
    switch (decoderType)
    {
	case decoderType::SCRecursive: {
		decoderPtr = new ScRecursiveDecoder(codePtr);
		break;
	}
	case decoderType::SC: {
		decoderPtr = new ScDecoder(codePtr);
	}
		break;
	case decoderType::SCFano: {
		double T = ExtractDouble(decoderParams, "T", "SCFano decoder");
		double delta = ExtractDouble(decoderParams, "delta", "SCFano decoder");
		decoderPtr = new ScFanoDecoder(codePtr, T, delta);
		return decoderPtr;
	}
	case decoderType::SCFlipFano: {
		double T = ExtractDouble(decoderParams, "T", "SCFanoList decoder");
		double delta = ExtractDouble(decoderParams, "delta", "SCFanoList decoder");
		double approximationSnr = ExtractDouble(decoderParams, "ApproximationSnr", "SCFanoList decoder");
		double L = ExtractInt(decoderParams, "L", "SCFanoList decoder");

		decoderPtr = new ScFlipFanoDecoder(codePtr, T, delta, approximationSnr, L);
		return decoderPtr;
	}
		break;
	case decoderType::SCFlip: {
		int T = ExtractInt(decoderParams, "T", "SCFlip decoder");
		decoderPtr = new ScFlipDecoder(codePtr, T);
		break;
	}
	case decoderType::SCList: {
		int L = ExtractInt(decoderParams, "L", "SCList decoder");
		decoderPtr = new ScListDecoder(codePtr, L);
	}
	break;
	case decoderType::SCFlipProg: {
		int level = ExtractInt(decoderParams, "level", "SCFlipProg decoder");
		int gammaLeft = ExtractInt(decoderParams, "gammaLeft", "SCFlipProg decoder");
		int gammaRight = ExtractInt(decoderParams, "gammaRight", "SCFlipProg decoder");
		std::string omegaArrString = ExtractString(decoderParams, "omegaArr", "SCFlipProg decoder", true);
		std::vector<double> omegaArr = OmegaArrStrToVector(omegaArrString);

		decoderPtr = new ScFlipProgDecoder(codePtr, level, gammaLeft, gammaRight, omegaArr);
	}
		break;
	case decoderType::SCListFlipStat: {
		int L = ExtractInt(decoderParams, "L", "SCListFlipStat decoder");
		double omega = ExtractDouble(decoderParams, "Omega", "SCListFlipStat decoder");
		decoderPtr = new ScListFlipStatDecoder(codePtr, L, omega);
	}
		break;
	case decoderType::SCListFlipOracleStat: {
		int L = ExtractInt(decoderParams, "L", "SCListFlipOracleStat decoder");
		decoderPtr = new ScListFlipOracleStatDecoder(codePtr, L);
	}
	break;
    default:
        break;
    }
    
    return decoderPtr;
}

BaseSimulator * BuildSimulator(
                               simulatorType simulationType,
                               std::unordered_map<std::string, std::string> simulationTypeParams,
							   PolarCode * codePtr,
							   Encoder * encoderPtr,
                               BaseDecoder * decoderPtr)
{
    BaseSimulator * simulator = NULL;
    
    switch (simulationType)
    {
    case simulatorType::MC: {
		bool isSigmaDependOnR = (bool)ExtractInt(simulationTypeParams, "isSigmaDependOnR", "MC Simulator");
        int maxTestsCount = ExtractInt(simulationTypeParams, "maxTestsCount", "MC simulator");
        int maxRejectionsCount = ExtractInt(simulationTypeParams, "maxRejectionsCount", "MC simulator");
		std::string additionalInfoFilename = ExtractString(simulationTypeParams, "additionalInfoFilename", "MC simulator", false);
		std::string channelStr = ExtractString(simulationTypeParams, "channel", "MC simmulator" , false);
		
		BaseChannel * channelPtr;
		if (channelStr == "BPSK-AWGN")
			channelPtr = new BpskAwgnChannel();
		else if (channelStr == "BSC")
			channelPtr = new BscChannel((double)codePtr->k() / codePtr->N());
		else
			throw UnknownChannelException("Unknown channel: " + channelStr + " in simulator params section");

        simulator = new MonteCarloSimulator(maxTestsCount, maxRejectionsCount, codePtr, encoderPtr, channelPtr, decoderPtr, isSigmaDependOnR);
    }
        break;
    default:
        break;
    }
    return simulator;
}

void LogIntoFile(std::string filename, std::string message, bool isRewrite=false, std::string stringPrefix="") {
	
	std::ofstream resultsFileStream;
	if (isRewrite)
		resultsFileStream.open(filename, std::fstream::out | std::fstream::trunc);
	else
		resultsFileStream.open(filename, std::fstream::out | std::fstream::app);

	if (!resultsFileStream.is_open()) {
		throw FileIsNotOpennedException("Cannot open file \"" + filename + "\".");
	}

	std::string sentense;
	std::istringstream splitStream(message);
	while (std::getline(splitStream, sentense, '\n'))
	{
		if (!stringPrefix.empty())
			resultsFileStream << stringPrefix;

		resultsFileStream << sentense << std::endl;
	}
	
	resultsFileStream.close();
}

bool TryLogIntoFile(std::string filename, std::string message, bool isRewrite=false, std::string stringPrefix = "") {
	try {
		LogIntoFile(filename, message, isRewrite, stringPrefix = "");
		return true;
	}
	catch (const std::exception&) {
		return false;
	}
}

bool IsFileExists(std::string filename) {
	FILE *file;
	fopen_s(&file, filename.c_str(), "r");
	if (file) {
		fclose(file);
		return true;
	}
	else {
		return false;
	}
}

void LogIntoConsole(std::string message) {
	std::cout << message;
}

void Simulate(std::string configFilename) {
	PolarCode * codePtr = NULL;
	Encoder * encoderPtr = NULL;
    BaseDecoder * decoderPtr = NULL;
    BaseSimulator * simulatorPtr = NULL;
	std::vector<SimulationParams> simulationParamsArray;
	
	try {
		LogIntoConsole("Simulation initializing is starting....\n");

		simulationParamsArray = ReadConfig(configFilename);
	}
	catch (const std::exception& err) {
		std::string message = "Error was ocurred:\n" + std::string(err.what()) + "\n";
		LogIntoConsole(message);
		return;
	}

	LogIntoConsole("\tConfiguration file has been read succesfully.\n");
	for (auto simulationParams : simulationParamsArray) {
		try {
			codePtr = BuildCode(simulationParams.codeParams);

			encoderPtr = new Encoder(codePtr);
			decoderPtr = BuildDecoder(simulationParams.decoder, simulationParams.decoderParams, codePtr);
			LogIntoConsole("\tDecoder has been built succesfully.\n");

			simulatorPtr = BuildSimulator(simulationParams.simulator, simulationParams.simulatorParams, codePtr, encoderPtr, decoderPtr);
			LogIntoConsole("\tSimulator has been built succesfully.\n");

			bool isResultFileRewrited = false;
			if(!IsFileExists(simulationParams.resultsFilename) || isResultFileRewrited)
				LogIntoFile(simulationParams.resultsFilename, SimulationIterationResults::GetHeader() + "\n", false);

			LogIntoConsole("Simulation has been started.\n\n");
			LogIntoFile(simulationParams.resultsFilename, simulationParams.ToString(), false, "# ");
			LogIntoConsole(simulationParams.ToString());

			for (size_t i = 0; i < simulationParams.snrArray.size(); i++)
			{
				LogIntoConsole("Iteration has been started. SNR: " + std::to_string(simulationParams.snrArray[i]) + "\n");

				auto result = simulatorPtr->Run(simulationParams.snrArray[i]);
				auto message = result.ToString() + "\n";

#ifdef DECODER_STAT
				auto decoderStat = decoderPtr->GetStatistic();
				LogIntoFile(simulationParams.additionalFilename, "Statistic for SNR: " + std::to_string(simulationParams.snrArray[i]), false);
				bool isFileRewriting = false;
				LogIntoFile(simulationParams.additionalFilename, decoderStat, isFileRewriting);
				decoderPtr->ClearStatistic();
#endif // DECODER_STAT

				LogIntoFile(simulationParams.resultsFilename, message, false);
				LogIntoConsole("Iteration has been ended with result:\n\t" + message);
			}
		}
		catch (const std::exception& err) {
			std::string message = "Error was ocurred:\n" + std::string(err.what()) + "\n";
			LogIntoConsole(message);
			if (!TryLogIntoFile(simulationParams.resultsFilename, message, false))
				LogIntoConsole("Cannot write message of error into file \"" + simulationParams.resultsFilename + "\".\n");
		}
    }

	delete simulatorPtr;
	delete decoderPtr;
	delete encoderPtr;
	delete codePtr;
}
