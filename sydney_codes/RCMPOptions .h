#ifndef RCMPOPTIONS_H
#define RCMPOPTIONS_H

struct RCMPOptions {

    // Directories
    static inline std::string const histogram_directory = "/home/sydney/S2232_analysis/Histograms/RCMP";
    static inline std::string const analysis_directory = "/media/sydney/GRINYER_28/S2232_analysis/AnalysisTrees/";
    static inline std::string const fragment_directory = "/media/sydney/GRINYER_28/S2232_analysis/FragmentTrees/";
    static inline std::string const calibration_directory = "/home/sydney/S2232_analysis/cal_files/";
    static inline std::string const calibration_directory_files = "/home/sydney/S2232_analysis/cal_files/files/";
    static inline std::string const default_calibration = "/home/sydney/S2232_analysis/cal_files/S2232Cal.cal/";
    static inline std::string const log_directory = "/home/sydney/S2232_analysis/Logs/EnergyCalLogs/ ";
    static inline std::string const decay_fit_directory = "/home/sydney/S2232_analysis/DecayFit/";
    static inline std::string const pile_up_fit_directory = "/home/sydney/S2232_analysis/PUFit/";
    static inline std::string const data_directory = "/media/sydney/RCMP_SG/S2232_midas_files/";

    // Histogram prefixes
    static inline std::string const rcmp_prefix = "RCMP";
    static inline std::string const analysis_prefix = "Analysis";
    static inline std::string const combined_prefix = "Combined";

    static inline int const first_run = 28183; //28183
    static inline int const last_run = 28351; //28351
    static inline int const matrix_dimension = 1024; // How do you get this number?
    static inline int const number_of_strips = 32;
    static inline int const number_of_detectors = 6;
    
};

#endif
