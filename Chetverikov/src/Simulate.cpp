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
#include "../include/ScDecoder.h"
#include "../include/ScFanoDecoder.h"
#include "../include/ScFlipDecoder.h"
#include "../include/BaseSimulator.h"
#include "../include/ConfigReading.h"
#include "../include/Exceptions.h"

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

std::vector<int> StrToVector(std::string crcPolyStr) {
	if (crcPolyStr.empty())
		return std::vector<int>();

	size_t deg = crcPolyStr.size();

	std::vector<int> result(deg + 1, 0);
	result[deg] = 1;

	for (size_t i = 0; i < deg; i++)
	{
		if (crcPolyStr[i] == '1')
			result[deg - i - 1] = 1;

		else if (crcPolyStr[i] != '0')
			throw CrcPolyException("Crc polynom has incorrect symbol: " + crcPolyStr);
	}

	return result;
}

PolarCode * BuildCode(std::unordered_map<std::string, std::string> codeParams) {
	PolarCode * codePtr;

	int m = ExtractInt(codeParams, "m", "PolarCode");
	int k = ExtractInt(codeParams, "k", "PolarCode");
	std::string crcPolyString = ExtractString(codeParams, "CRC", "PolarCode", false);
	std::vector<int> crcPoly = StrToVector(crcPolyString);
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
	case decoderType::SC: {
		decoderPtr = new ScDecoder(codePtr);
	}
		break;
	case decoderType::SCFano: {
		double T = ExtractDouble(decoderParams, "T", "SCFano decoder");
		double delta = ExtractDouble(decoderParams, "delta", "SCFano decoder");
		decoderPtr = new ScFanoDecoder(codePtr, T, delta);
	}
	case decoderType::SCFlip: {
		int T = ExtractInt(decoderParams, "T", "SCFano decoder");
		decoderPtr = new ScFlipDecoder(codePtr, T);
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
        int maxTestsCount = ExtractInt(simulationTypeParams, "maxTestsCount", "MC simulator");
        int maxRejectionsCount = ExtractInt(simulationTypeParams, "maxRejectionsCount", "MC simulator");
            
        simulator = new MonteCarloSimulator(maxTestsCount, maxRejectionsCount, codePtr, encoderPtr, decoderPtr);
    }
        break;
    default:
        break;
    }
    return simulator;
}

void LogIntoFile(std::string filename, std::string message, std::string stringPrefix="") {
	
	std::ofstream resultsFileStream;
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

bool TryLogIntoFile(std::string filename, std::string message, std::string stringPrefix = "") {
	try {
		LogIntoFile(filename, message, stringPrefix = "");
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

			if(!IsFileExists(simulationParams.resultsFilename))
				LogIntoFile(simulationParams.resultsFilename, SimulationIterationResults::GetHeader() + "\n");

			LogIntoConsole("Simulation has been started.\n\n");
			LogIntoFile(simulationParams.resultsFilename, simulationParams.ToString(), "# ");
			LogIntoConsole(simulationParams.ToString());


			for (size_t i = 0; i < simulationParams.snrArray.size(); i++)
			{
				LogIntoConsole("Iteration has been started. SNR: " + std::to_string(simulationParams.snrArray[i]) + "\n");

				auto result = simulatorPtr->Run(simulationParams.snrArray[i]);
				auto message = result.ToString() + "\n";

				LogIntoFile(simulationParams.resultsFilename, message);
				LogIntoConsole("Iteration has been ended with result:\n\t" + message);
			}
		}
		catch (const std::exception& err) {
			std::string message = "Error was ocurred:\n" + std::string(err.what()) + "\n";
			LogIntoConsole(message);
			if (!TryLogIntoFile(simulationParams.resultsFilename, message))
				LogIntoConsole("Cannot write message of error into file \"" + simulationParams.resultsFilename + "\".\n");
		}
    }

	delete simulatorPtr;
	delete decoderPtr;
	delete codePtr;
}
