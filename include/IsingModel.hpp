#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

#include "utils.hpp"



// GRID INIT FLAGS
enum eInitFlags {
    INIT_T_0,
    INIT_T_INF
};

// SIMULATION TYPES
enum eSimTypes {
    SINGLE_SIM,
    MULTI_SIM,
    RAYLIB_SIM
};


// SIMULATION INFO
typedef struct {
    vector<double> temps;

    // calculated during single_sims
    vector<double> mags;
    vector<double> energies;
    
    // calculated during multi_sims
    vector<double> susceptibilities;
    vector<double> specific_heats;
} IsingModelInfo;



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

    TFile *hfile;
    TNtuple *data;

    IsingModelInfo ismdl_info;
public:
    
    /* ************* INITIALIZATION FUNCTIONS *********** */
    IsingModel();
    ~IsingModel();


    inline void set_temp(double _temp) {beta = 1.0/_temp;}
    void init_grid(int flag);


    int calc_total_energy();
    double calc_total_energy_per_site();
    int calc_total_mag();
    double calc_total_mag_per_site();
    int calc_delta_U(int i, int j); // calculates energy of set of dipole at (i,j) and xes/xer 4 neighbors

    
    void update_grid();
    int64_t simulate(int N_iter);   // runs one simulation at one temperature
    int64_t multi_simulate(         // runs multiple simulations in a range of temperatures
        double T0,
        double Tf,
        double dT,
        double N_iter_equib,
        double N_iter_avg
    );
    void show_simulation(
        int win_w, int win_h,
        string win_title,
        int refresh_rate
    );


    void run();
    
};




#endif // ISINGMODEL_HPP