#include "IsingModel.hpp"
#include "Constants.hpp"


IsingModel::IsingModel() : grid(nullptr), hfile(nullptr), data(nullptr) {

    // grab all data from the configuration file
    ifstream f;
    f.open(CONFIG_FILE_PATH);

    // store all data in a map
    string line;
    while (getline(f, line)) {
        if ((line.size() < 1) || (line[0] == '#')) continue;

        size_t equals_pos = line.find('=');
        string key = line.substr(0, equals_pos);
        string value = line.substr(equals_pos+1, line.size() - equals_pos);

        init_data[key] = value;
    }
    f.close();

    // now we need to set the settings that are universal for any simulation
    // this is: lattice size, initial temperature


    // first find the simulation type
    sim_type = stod(init_data["SIM_TYPE"]);
    lattice_size = stoi(init_data["LATTICE_SIZE"]);

    // initilize the lattice with the chosen size
    grid = new int*[lattice_size];
    for (int i=0; i<lattice_size; i++) {
        grid[i] = new int[lattice_size];
        for (int j=0; j<lattice_size; j++) {
            grid[i][j] = 0;
        }
    }

    // set and store the initial temp
    // needs to be stored for raylib simulation, since I want to be able to press R and restart the simulation
    grid_init_flag = stoi(init_data["INITIAL_TEMP"]);
    InitGrid(grid_init_flag);

    cout << "Model Initialized\n";
}

IsingModel::~IsingModel() {
    // free memory
    for (int i=0; i<lattice_size; i++) {
        delete[] grid[i];
    }
    delete[] grid;

    if (hfile) delete hfile;
}





double IsingModel::CalcEnergy() {
    if (total_energy != 0) return total_energy; // it will have already been calculated

    // otherwise, we need to calculate it! this only is reached in a call to init_grid()
    int energy = 0;
    for (int i=0; i<lattice_size; i++) {
        for (int j=0; j<lattice_size; j++) {
            int i_minus1 = (i-1 < 0) ? lattice_size-1 : i-1;
            int j_minus1 = (j-1 < 0) ? lattice_size-1 : j-1;
            int i_plus1 = (i+1 > lattice_size-1) ? 0 : i+1;
            int j_plus1 = (j+1 > lattice_size-1) ? 0 : j+1;
            
            int energy_neighbors = -(grid[i_plus1][j] + grid[i_minus1][j] + grid[i][j_plus1] + grid[i][j_minus1]);
            int energy_current = grid[i][j]*energy_neighbors;

            energy += energy_current;
        }
    }
    energy = (int)energy/2; // the above method double counts; this fixes that
    total_energy = energy; // set the class total_energy so that it can be returned for other functions

    return (double)energy/(lattice_size*lattice_size);;
}



double IsingModel::CalcMagnetization() {
    int mag = 0;
    for (int i=0; i<lattice_size; i++) {
        for (int j=0; j<lattice_size; j++) {
            mag += grid[i][j];
        }
    }
    return (double)mag/(lattice_size*lattice_size);
}


int IsingModel::CalcDeltaE(int i, int j) {
    int i_minus1 = (i-1 < 0) ? lattice_size-1 : i-1;
    int j_minus1 = (j-1 < 0) ? lattice_size-1 : j-1;
    int i_plus1 = (i+1 > lattice_size-1) ? 0 : i+1;
    int j_plus1 = (j+1 > lattice_size-1) ? 0 : j+1;
    
    int energy_neighbors = -(grid[i_plus1][j] + grid[i_minus1][j] + grid[i][j_plus1] + grid[i][j_minus1]);
    int energy_current = grid[i][j]*energy_neighbors;
    int energy_proposed = -energy_current;
    int delta_U = energy_proposed - energy_current;

    return delta_U;
}



void IsingModel::InitGrid(int flag) {
    // when we init or reset the grid, we want to also reset the total energy
    // so that it can be calculated again.
    total_energy = 0; 

    if (flag == INIT_T_0) {
        for (int i=0; i<lattice_size; i++) {
            for (int j=0; j<lattice_size; j++) {
                grid[i][j] = 1;
            }
        }
    } else if (flag == INIT_T_INF) {
        for (int i=0; i<lattice_size; i++) {
            for (int j=0; j<lattice_size; j++) {
                double r = rand_double();
                grid[i][j] = (r < 0.5) ? 1 : -1;
            }
        }
    }
}

void IsingModel::UpdateGrid() {

    for (size_t n=0; n<(lattice_size*lattice_size); n++) {
        // pick random lattice point
        int i = rand()%lattice_size, j = rand()%lattice_size;

        // compute potential change in energy
        int delta_U = CalcDeltaE(i, j);
        

        if (delta_U <= 0) {         // if change results in lower (or zero) energy, flip automatically
            grid[i][j] *= -1;
            total_energy += delta_U;
        } else {                    // otherwise, we use the exponential Boltzmann factor as a probability to accept anyway
            double fact = exp(-beta*delta_U);
            double r = rand_double();
            if (r < fact) {
                grid[i][j] *= -1;
                total_energy += delta_U;
            }
        }
    }
}


int64_t IsingModel::Simulate(int N_iter) {
    auto start = chrono::high_resolution_clock::now();

    for (int i=0; i<N_iter; i++) {
        UpdateGrid();
        double energy = CalcEnergy();
        double mag = CalcMagnetization();
        data->Fill((float)i, energy, mag);
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);

    return duration.count();
}

int64_t IsingModel::MultiSimulate(double T0, double Tf, double dT, double N_iter_equib, double N_iter_avg) {

    // start timer
    auto start = chrono::high_resolution_clock::now();

    int N_iter = (int)((Tf-T0)/dT) + 2; // +1 from double->int, +1 from exclusive upper bound on loop
    for (int n=0; n<N_iter; n++) {

        // calculate and set temperature for this iteration
        double temp = T0 + n*dT;
        SetTemp(temp);

        // first loop gets system in equilibrium
        for (int i=0; i<N_iter_equib; i++) {
            UpdateGrid();
        }

        // second loop is where values/averages are computed.
        double avg_mag = 0.0;
        double avg_energy = 0.0;
        double mag_squared = 0.0;
        double energy_squared = 0.0;
        for (int i=0; i<N_iter_avg; i++) {
            double energy = CalcEnergy();
            double mag = CalcMagnetization();
            avg_mag += abs(mag); // absolute value of the mean magnetization for a given run
            avg_energy += energy;
            mag_squared += mag*mag;
            energy_squared += energy*energy;

            UpdateGrid();
        }
        // divide the quantities by the number of samples that were taken to get means
        avg_mag /= N_iter_avg;
        avg_energy /= N_iter_avg;
        mag_squared /= N_iter_avg;
        energy_squared /= N_iter_avg;


        // add all of the values to vectors to plot

        double suscept = beta*lattice_size*(mag_squared - avg_mag*avg_mag);
        double heat = beta*beta*(energy_squared - avg_energy*avg_energy)/lattice_size;

        data->Fill(temp, avg_energy, avg_mag, suscept, heat);

        ResetEnergy(); // reset the total energy for the next simulation
        cout << "T=" << temp << " sim complete." << "\n"; // good way to know that it is working!
    }

    auto end = chrono::high_resolution_clock::now();
    cout << "All sims complete\n";

    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start); // prob should return seconds...
    return duration.count();
}

void IsingModel::ShowSimulation(
    int win_w, int win_h,
    string win_title,
    int refresh_rate
) {

    int n = lattice_size; // for easier typing
    InitWindow(win_w, win_h, win_title.c_str());
    SetTargetFPS(refresh_rate);

    while (!WindowShouldClose()) {
        BeginDrawing();

        // want to be able to restart it if we press "r"
        if (IsKeyPressed(KEY_R)) {
            InitGrid(grid_init_flag);
        }

        ClearBackground(BLACK);
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                Color color = (grid[i][j] > 0) ? WHITE : BLACK;
                DrawRectangle(j*(win_w/n), i*(win_h/n), (win_w/n), (win_h/n), color);
            }
        }
        UpdateGrid();

        EndDrawing();
    }

    CloseWindow();
}




void IsingModel::Run() {
    // each part of the if(else) statement just grabs the values specific for that run,
    // then executes the corresponding function

    if (sim_type == SINGLE_SIM) {
        double equib_temp = stod(init_data["EQUIB_TEMP"]);
        SetTemp(equib_temp);
    
        // if we have equilibration time override or not
        int equib_override = stoi(init_data["OVERRIDE_EQUIB_TIME"]);
        int N_iter = (equib_override) ? stoi(init_data["EQUIB_TIME"]) : NUM_EQUIB_ITER;

        // create the data ntuple and file handle for a single sim
        // test if the file already exists, if so, basically just delete/clear it and recreate it
        const TString file_path = DATA_FILE_PATH + SINGLE_SIM_FILE_NAME;
        hfile = (TFile *)gROOT->FindObject(file_path); if (hfile) hfile->Close();
        hfile = new TFile(file_path, "RECREATE", "Single simuation data");
        // create the ntuple to contain the iterations, energies, and magnetizations
        data = new TNtuple("data", "Single simulation data", "N:energy:mag");

        // perform the simulation
        cout << "Running single-simulation...\n";
        int64_t ms = Simulate(N_iter);
        cout << "\nDone! Simulation took " << ms << " milliseconds\n";

        // write the ntuple data to the file
        data->Write();
    } else if (sim_type == MULTI_SIM) {
        double t0 = stod(init_data["T0"]);
        double tf = stod(init_data["TF"]);
        double dt = stod(init_data["DT"]);

        // create the file handle for the multi simulation data
        const TString file_path = DATA_FILE_PATH + MULTI_SIM_FILE_NAME;
        hfile = (TFile *)gROOT->FindObject(file_path); if (hfile) hfile->Close();
        hfile = new TFile(file_path, "RECREATE", "Multi simulation data");
        // now we create the ntuple
        data = new TNtuple("data", "Multi simulation data", "T:energy:mag:suscept:heat");

        cout << "Running multi-simulation...\n";
        int64_t ms = MultiSimulate(t0, tf, dt, NUM_EQUIB_ITER, NUM_AVG_ITER);
        cout << "\nDone! All simulations took " << ms << " milliseconds, or " << (double)ms/1000.0 << " seconds.\n";

        // after the simulation finishes, we write it to the file
        data->Write();
    } else if (sim_type == RAYLIB_SIM) {

        int win_w = stoi(init_data["WIN_W"]);
        int win_h = stoi(init_data["WIN_H"]);
        string win_title = init_data["WIN_TITLE"];
        int refresh_rate = stoi(init_data["REFRESH_RATE"]);

        int equib_temp = stod(init_data["EQUIB_TEMP"]);
        SetTemp(equib_temp);

        ShowSimulation(win_w, win_h, win_title, refresh_rate);
    }
}