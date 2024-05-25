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

/* ********* ISING MODEL ************* */
class IsingModel {
private:
    // lattice info
    int lattice_size;
    int **grid;
    int grid_init_flag;
    double beta;
    int sim_type; // the type of simulation currently running

    int total_energy; // current total energy in the lattice
    map<string, string> init_data; // all the data from run.config

    // file handle and data ntuple for data storage
    TFile *hfile;
    TNtuple *data;
public:
    
    /* ************* INITIALIZATION FUNCTIONS *********** */
    IsingModel();
    ~IsingModel();


    inline void SetTemp(double _temp) {beta = 1.0/_temp;}
    void InitGrid(int flag);
    void inline ResetEnergy() {total_energy = 0;} // for clarity inside of functions

    double CalcEnergy();
    double CalcMagnetization();
    int CalcDeltaE(int i, int j); // calculates energy of set of dipole at (i,j) and xes/xer 4 neighbors

    
    void UpdateGrid();
    int64_t Simulate(int N_iter);   // runs one simulation at one temperature
    int64_t MultiSimulate(         // runs multiple simulations in a range of temperatures
        double T0,
        double Tf,
        double dT,
        double N_iter_equib,
        double N_iter_avg
    );
    void ShowSimulation(
        int win_w, int win_h,
        string win_title,
        int refresh_rate
    );


    void Run();
    
};




#endif // ISINGMODEL_HPP