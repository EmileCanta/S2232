#ifndef FINDHISTOGRAMFILE_H
#define FINDHISTOGRAMFILE_H

#include <filesystem>

std::string findHistogramFile(int t_run_number, std::string t_histogram_directory, std::string t_histogram_prefix) {

    std::string file_prefix = t_histogram_prefix+std::to_string(t_run_number);
    std::string directory_file;
    std::string file_name;
    int file_hits = 0;
	for (const auto &entry : std::filesystem::directory_iterator(t_histogram_directory)) {
        directory_file = entry.path();
        if (directory_file.find(file_prefix) != std::string::npos) {
            if (file_hits == 0) {
				file_name = entry.path();
			}
			file_hits++;
        }
    }
	if (file_hits == 0) {
		std::cerr << "[WARNING] No histogram found for run " << t_run_number << "!" << std::endl;
        return "";
	}
	else if (file_hits > 1) {
		std::cerr << "[WARNING] More than one histogram found for run " << 
		t_run_number << "; using first histogram found." << std::endl;
	}

    return file_name;
}

#endif