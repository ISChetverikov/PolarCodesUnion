
#include <unordered_map>
#include "../include/SimulationParameters.h"
#include "../include/ConfigReading.h"
#include "../lib/pugixml-1.10/src/pugixml.hpp"
#include "../lib/pugixml-1.10/src/pugiconfig.hpp"
#include "../include/Exceptions.h"
#include "../include/SimulationType.h"
#include "../include/DecoderType.h"

#define ROOT_SECTION "Config"
#define		SIMULATION_SECTION "Simulation"
#define			RESULT_FILE_RECORD "ResultFile"
#define			ADDITIONAL_FILE_RECORD "AdditionalFile"
#define			SIMULATOR_SECTION "Simulator"
#define			CODE_SECTION "Code"
#define			DECODER_SECTION "Decoder"
#define			SNR_RANGE_SECTION "SnrRange"
#define				SNR_RANGE_START_ATTRIBUTE "Start"
#define				SNR_RANGE_STOP_ATTRIBUTE "Stop"
#define				SNR_RANGE_STEP_ATTRIBUTE "Step"
#define			SNR_ARRAY_SECTION "SnrArray"
#define				SNR_RECORD "Snr"

#define NAME_ATTRIBUTE "name"
#define TYPE_ATTRIBUTE "type"
#define PARAM_RECORD "Param"

// if not exist then exception
pugi::xml_node GetChildNode(pugi::xml_node node, std::string childName) {
	auto child = node.child(childName.c_str());
	if (child == NULL)
		throw ConfigParseException("Parse config error: can not find \"" + childName + "\" section");

	return child;
}

// if not exists then returns NULL
pugi::xml_node TryGetChildNode(pugi::xml_node node, std::string childName) {
	auto child = node.child(childName.c_str());
	
	return child;
}

pugi::xml_attribute GetAttribute(pugi::xml_node node, std::string name) {
	auto attr = node.attribute(name.c_str());
	if (attr == NULL)
		throw ConfigParseException("Parse config error: can not find \"" + name + "\" atribute");

	return attr;
}

std::unordered_map <std::string, std::string > ReadParams(pugi::xml_node node) {
	std::unordered_map <std::string, std::string > params;
	for (pugi::xml_node param_node : node.children(PARAM_RECORD))
	{
		params[GetAttribute(param_node, NAME_ATTRIBUTE).value()] = param_node.child_value();
	}

	return params;
}

SimulationParams ReadSimulationSection(pugi::xml_node simulation_node) {
	SimulationParams params;

	params.resultsFilename = GetChildNode(simulation_node, RESULT_FILE_RECORD).child_value();
	params.additionalFilename = GetChildNode(simulation_node, ADDITIONAL_FILE_RECORD).child_value();

	pugi::xml_node simulator_node = GetChildNode(simulation_node, SIMULATOR_SECTION);

	std::string simulatorStr = GetAttribute(simulator_node, TYPE_ATTRIBUTE).value();
	simulatorType simulator = simulatorFromString(simulatorStr);
	if (simulator == simulatorType::UnknownSimulation)
		throw ConfigParseException("Parse config error: unsupported value of simulation type: " + simulatorStr);

	params.simulator = simulator;
	params.simulatorParams = ReadParams(simulator_node);

	pugi::xml_node code_node = GetChildNode(simulation_node, CODE_SECTION);

	std::string codeStr = GetAttribute(code_node, TYPE_ATTRIBUTE).value();
	params.code = codeStr;
	params.codeParams = ReadParams(code_node);

	pugi::xml_node decoder_node = GetChildNode(simulation_node, DECODER_SECTION);

	std::string decoderStr = GetAttribute(decoder_node, TYPE_ATTRIBUTE).value();
	decoderType decoder = decoderTypeFromString(decoderStr);
	if (decoder == decoderType::UnknownDecoder)
		throw ConfigParseException("Parse config error: unsupported value of decoder type: " + decoderStr);

	params.decoder = decoderTypeFromString(decoderStr);
	params.decoderParams = ReadParams(decoder_node);

	pugi::xml_node snrRange_node = TryGetChildNode(simulation_node, SNR_RANGE_SECTION);
	bool isSnrPresent = false;
	if (snrRange_node != NULL) {
		for (double i = std::stod(GetAttribute(snrRange_node, "start").value());
			i < std::stod(GetAttribute(snrRange_node, "stop").value());
			i += std::stod(GetAttribute(snrRange_node, "step").value()))
			params.snrArray.push_back(i);
		isSnrPresent = true;
	}

	pugi::xml_node snrArray_node = TryGetChildNode(simulation_node, SNR_ARRAY_SECTION);
	if (snrArray_node != NULL) {
		for (pugi::xml_node snr_node : snrArray_node.children(SNR_RECORD))
		{
			params.snrArray.push_back(std::stod(snr_node.child_value()));
		}
		isSnrPresent = true;
	}

	if (!isSnrPresent)
		throw ConfigParseException("Parse config error: either \"SnrArray\" and \"SnrRange\" sections are absent.");

	return params;
}

std::vector<SimulationParams> ReadConfig(std::string configFilename) {
	std::vector<SimulationParams> paramsArray;

	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(configFilename.c_str());

	if (!result)
	{
		throw ConfigParseException("Error during open config xml file:" + std::string(result.description()));
	}

	pugi::xml_node config_node = GetChildNode(doc, ROOT_SECTION);
	for (pugi::xml_node simulation_node : config_node.children(SIMULATION_SECTION))
	{
		paramsArray.push_back(ReadSimulationSection(simulation_node));
	}

	if (paramsArray.empty())
		throw ConfigParseException("Parse config error: any \"Simulation\" section is absent.");

	return paramsArray;
}