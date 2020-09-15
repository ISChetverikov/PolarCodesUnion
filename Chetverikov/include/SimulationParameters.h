#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <iostream>
#include <vector>
#include "DecoderType.h"
#include "SimulationType.h"

struct SimulationParams {
    simulatorType simulator;
    std::unordered_map < std::string, std::string > simulatorParams;
    
	decoderType decoder;
	std::unordered_map < std::string, std::string > decoderParams;
	
	// Others code types are not forseen yet
	std::string code = "Polar";
	std::unordered_map < std::string, std::string > codeParams;

	std::string resultsFilename;
	std::string additionalFilename;
	std::vector<double> snrArray;

	std::string ToString() {
		std::stringstream ss;

		ss << std::string(60, '#') << "\n";
		ss << "Prameters of simulation run\n\n";

		ss << "ResultsFilename: " + resultsFilename + "\n\n";
		ss << "AdditionalFilename: " + additionalFilename + "\n\n";
		
		ss << "Simulation type: " + simulationTypeToString(simulator) + "\n";
		ss << "Simulation parameters:\n";
		for (auto simulationParam : simulatorParams)
		{
			ss << "\t" << simulationParam.first + ": " + simulationParam.second + "\n";
		}

		ss << "Decoder type: " + decoderTypeToString(decoder) + "\n";
		ss << "Decoder params:\n";
		for (auto decoderParam : decoderParams)
		{
			ss << "\t" << decoderParam.first + ": " + decoderParam.second + "\n";
		}

		ss << "Code type: " + code + "\n";
		ss << "Code params:\n";
		for (auto codeParam : codeParams)
		{
			ss << "\t" << codeParam.first + ": " + codeParam.second + "\n";
		}

		ss << "\n";
		ss << "SNR: ";
		for (size_t i = 0; i < snrArray.size(); i++)
		{
			ss <<  snrArray[i] << ", ";
		}
		ss << "\n";

		ss << std::string(60, '#') << "\n";

		return ss.str();
	}
};

struct SimulationIterationResults {
	double snr;
	double ebn0;
	double sigma;

	double fer;

	int rejectionsCount;
	int testsCount;
	std::chrono::milliseconds elapsedTime;

	static std::string GetHeader() {
		return "SNR, EbN0, sigma, FER, rejectionsCount, testsCount, time(ms)";
	}

	std::string ToString() {
		std::stringstream ss;
		
		ss << snr << ", " << ebn0 << ", " << sigma << ", " << fer << ", " << rejectionsCount
			<< ", " << testsCount << ", " << elapsedTime.count();
		
		return ss.str();
	}
};