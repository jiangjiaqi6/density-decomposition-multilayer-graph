#ifndef UTILITY_H_
#define UTILITY_H_

#include "Defines.h"


class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<ui> sequence;
		if(number == 0) sequence.push_back(0);
		while(number > 0) {
			sequence.push_back(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(unsigned int i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}

	static size_t getPeakRSS() {
        std::ifstream file("/proc/self/status");
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open /proc/self/status" << std::endl;
            return 0;
        }
		size_t peakRSS = 0;
		for (std::string line; std::getline(file, line);) {
			if (line.rfind("VmPeak:", 0) == 0) {
				peakRSS = std::stod(line.substr(7));
				std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
			}
			// if (line.rfind("VmSize:", 0) == 0) {
			// std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
			// }
			if (line.rfind("VmHWM:", 0) == 0) {
				std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
			}
			// if (line.rfind("VmRSS:", 0) == 0) {
			// std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
			// }
		}

        return 0; // VmHWM is reported in kB, convert to bytes
    }
	static size_t getSysPeakRSS() {
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) == 0) {
			return usage.ru_maxrss * 1024; // 返回的值是KB，乘以1024转换成字节
		} else {
			std::cerr << "Error: Unable to get resource usage" << std::endl;
			return 0;
		}
	}
};

#endif /* UTILITY_H_ */