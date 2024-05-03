#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

#include "utils.hpp"



// GRID INIT FLAGS
#define INIT_T_0 0
#define INIT_T_INF 1

// FILE PATHS
#define CONFIG_FILE_PATH "res/run.config"
#define GRAPH_OUTPUT_PATH "res/graphs"
#define SCREENSHOT_OUTPUT_PATH1 "/home/champson/Documents/School/Spring24/Comp/Project/src/res/screenshots"
#define SCREENSHOT_OUTPUT_PATH "./res/screenshots"

// SIMULATION TYPES
#define SINGLE_SIM 0
#define MULTI_SIM 1
#define RAYLIB_SIM 2

// DEFAULT SIMULATION VALUES
#define NUM_EQUIB_ITER 3000
#define NUM_AVG_ITER 15000


// OTHER DEFAULT VALS FOR CALCULATIONS
#define CORRELATION_TIME_RANGE 200
#define CORRELATION_TIME 200
#define T_EXACT 2.269


// SIMULATION INFO
typedef struct {
    vector<double> temps;

    // calculated during single_sims
    vector<double> mags;
    vector<double> energies;
    vector<double> autocorrelation;
    double auto_correlation_slope;
    double auto_correlation_y_int;
    
    // calculated during multi_sims
    vector<double> susceptibilities;
    vector<double> specific_heats;
    vector<double> correlation_time;
    vector<double> intermediate_mags; // used for autocorrelation calculation
} IsingModelInfo;



// RAYLIB INFO
#define MAX_SCREENSHOTS 10
typedef struct {
    // information needed for raylib
    int win_h;
    int win_w;
    string win_title;
    int refresh_rate;

    // information needed for showing the simulation
    int grid_init_flag;
    double equib_temp;

    // list of N's at which to take screenshots
    int screenshots[10];
} RaylibInfo;




/* ********* ISING MODEL ************* */
class IsingModel {
private:
    // lattice info
    int lattice_size;
    int **grid;
    int grid_init_flag;
    double beta;
    
    int sim_type; // the type of simulation currently running

    int total_energy; // the stored value of the total energy

    map<string, string> init_data; // all the data from run.config

    // the two structures containing further info that can be separated
    RaylibInfo raylib_info;
    IsingModelInfo ismdl_info;

public:
    
    /* ************* INITIALIZATION FUNCTIONS *********** */
    IsingModel();
    ~IsingModel();

    /**
     * Set's the equilibrium temperature. NOTE THAT THIS IS THE TEMPERATURE VALUE, NOT THE BETA VALUE
    */
    void set_temp(double _temp);

    /**
     * Initializes the grid to either T=infinity, which will randomly distribute spin up/down dipoles throughout,
     * or T=0, which makes every dipole spin up
    */
    void init_grid(int flag);




    /* ************ CALCULATION FUNTIONS *************** */

    /**
     * Essentially, this functions is called once at the beginning of a simulation,
     * after which it will set a class variable with the energy it calculates.
     * When it is called again, such as when averaging data during a simulation,
     * it will not recalculate anything, instead it will just return the class variable.
     * 
     * This probably could have been implemented a little bit better, but it works.
    */
    int calc_total_energy();

    /**
     * Returns the energy divided by the number of dipoles (N^2) i.e. the energy per lattice site
    */
    double calc_total_energy_per_site();

    /**
     * Calculates the total magnetization of the system by simply summing over all the spins.
    */
    int calc_total_mag();

    /**
     * Just like with `calc_total_energy_per_site()`, this returns the magnetization per lattice site.
    */
    double calc_total_mag_per_site();

    /**
     * This calculates autocorrelation function for the magnetization, and returns the correlation time.
     * It takes in the number of iterations used to reach equilibrium, 
     * as it starts sampling magnetization values after that.
     * The autocorrelation values are stored in the `ismdl_info` structure.
    */
    double calc_autocorrelation();

    /**
     * Given a lattice site (i,j), this calculates the potential change in energy
     * due to flipping that dipole.
    */
    int calc_delta_U(int i, int j); // calculates energy of set of dipole at (i,j) and xes/xer 4 neighbors

    


    /* *********** SIMULATION FUNCTIONS ************ */

    /**
     * Updates the grid once per lattice site (statistically), so N*N times.
    */
    void update_grid();

    /**
     * This runs a full, single simulation. It starts the model at the input temperature, either zero or infinity,
     * then lets it run either the defualt number of times, which is enough plus a little extra for a 100x100
     * lattice to equilibrate, or a different value based on settings in `run.config`.
     * 
     * It will compute the magnetization, the energy, and the autocorrelation after N iterations per lattice site.
     * 
     * @returns The time, in milliseconds, that it took to run the simulation
    */
    int64_t simulate(int N_iter);   // runs one simulation at one temperature

    /**
     * This runs multiple single simulations at a variety of temperatures designated in the `run.config` file.
     * 
     * This will compute the magnetization, specific heat, susceptibility, and autocorrelation per temperature.
     * 
     * @returns The time, in milliseconds, that it took to run the simulation.
     * This could certainly be converted into seconds; I just wanted to keep consistency.
    */
    int64_t multi_simulate(         // runs multiple simulations in a range of temperatures
        double T0,
        double Tf,
        double dT,
        double N_iter_equib,
        double N_iter_avg
    );


    /**
     * This generates all of the plots for the simulation. 
     * The paths that these are all saved to are relative, meaning they should work on any system,
     * but if they need to be changed, the top of this file contains defines for those paths
    */
    void plots();       // Creates and saves plots.


    /* ********** RAYLIB STUFF ************* */

    /**
     * Draws the grid, and colors spin-up dipoles as white, and spin-down dipoles as black.
    */
    void draw_grid();

    /**
     * The main render loop. Also takes and saves the requested screenshots.
     * The path that they are saved to is also relative, so just like with the plotting, it should work,
     * but the define for the file paths is at the top of this file if they need to be changed.
    */
    void show_simulation();


    /**
     * Big red button! Click at your own risk... 
     * Well, you kinda need to press this, if you don't then nothing will happen :-)
    */
    void run();
    
};




#endif // ISINGMODEL_HPP