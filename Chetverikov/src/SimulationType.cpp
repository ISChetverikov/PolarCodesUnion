#include "../include/SimulationType.h"

simulatorType simulatorFromString(std::string str){
	std::unordered_map<std::string, simulatorType> simulationTypeResolver = {
		{"MC", MC}
	};

	if (simulationTypeResolver.count(str) > 0)
		return simulationTypeResolver[str];
	
	return UnknownSimulation;
}

std::string simulationTypeToString(simulatorType type) {
	std::unordered_map<simulatorType, std::string> simulationTypeStringResolver = {
		{MC, "MC"},
		{UnknownSimulation, "UnknownSimulation"}
	};

	return simulationTypeStringResolver[type];
}

