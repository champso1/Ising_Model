#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <string>

// FILE PATHS
static const std::string CONFIG_FILE_PATH = "./res/run.config";
static const std::string DATA_FILE_PATH = "./res/plot_data/";
static const std::string SINGLE_SIM_FILE_NAME = "single_sim.root";
static const std::string MULTI_SIM_FILE_NAME = "multi_sim.root";

// DEFAULT SIMULATION VALUES
#define NUM_EQUIB_ITER 5000
#define NUM_AVG_ITER 15000

// OTHER DEFAULT VALS FOR CALCULATIONS
#define T_EXACT 2.269


#endif // CONSTANTS_HPP