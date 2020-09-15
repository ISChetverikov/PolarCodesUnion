#include <iostream>
#include <vector>
#include "../include/Simulate.h"

void PrintUsage() {
#ifdef __linux__ 
	std::cout << "Usage: LDPC config_path" << std::endl;
#elif _WIN32
	std::cout << "Usage: LDPC.exe config_path" << std::endl;
#else

#endif
}

int main(int argc, char* argv[]) {

	if (argc != 2) {
		PrintUsage();
		return 1;
	}
		
	Simulate(argv[1]);


}