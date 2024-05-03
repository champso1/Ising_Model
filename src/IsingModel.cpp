#include "IsingModel.hpp"


IsingModel::IsingModel() : grid(nullptr) {

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
    init_grid(grid_init_flag);

    cout << "Model Initialized\n";
}

IsingModel::~IsingModel() {
    // free memory
    for (int i=0; i<lattice_size; i++) {
        delete[] grid[i];
    }
    delete[] grid;
}





int IsingModel::calc_total_energy() {
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
    return energy;
}


double IsingModel::calc_total_energy_per_site() {
    return (double)calc_total_energy()/(lattice_size*lattice_size);
}


int IsingModel::calc_total_mag() {
    int mag = 0;
    for (int i=0; i<lattice_size; i++) {
        for (int j=0; j<lattice_size; j++) {
            mag += grid[i][j];
        }
    }
    return mag;
}



double IsingModel::calc_total_mag_per_site() {
    return (double)calc_total_mag()/(lattice_size*lattice_size);
}


int IsingModel::calc_delta_U(int i, int j) {
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



void IsingModel::init_grid(int flag) {
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



void IsingModel::set_temp(double _temp) {
    beta = 1.0/_temp;
}



double IsingModel::calc_autocorrelation() {
    // we grab the last 500 values from the magnetization to ensure that the model has equilibrated
    // we choose a good chunk of time, but only a few hundred, which is a few times the expected correlation time
    vector<double> autocorr_mags;
    if (sim_type == SINGLE_SIM) {
        for (int i=0; i<(CORRELATION_TIME_RANGE*4); i++) {
            autocorr_mags.push_back(ismdl_info.mags.at(i));
        }
    } else if (sim_type == MULTI_SIM) {
        ismdl_info.autocorrelation.clear(); // in a multisim, this may have been used before; we want to make sure it's empty before continuing
        for (int i=0; i<(CORRELATION_TIME_RANGE*4); i++) {
            autocorr_mags.push_back(ismdl_info.intermediate_mags.at(i));
        }
    }
    // now that we have grabbed the data from intermediate_mags, we can clear it
    ismdl_info.intermediate_mags.clear();


    //for_each(autocorr_mags.begin(), autocorr_mags.end(), [](double x){cout << x << " ";});
    int t_max = autocorr_mags.size();

    // we simply calculate using the formula.
    for (int t=0; t<t_max; t++) {
        double f0 = (double)1.0/(t_max-t);
        double factor1 = 0.0;
        double factor2 = 0.0;
        double factor3 = 0.0;
        for (int tt=0; tt<(t_max-t); tt++) {
            factor1 += autocorr_mags.at(tt)*autocorr_mags.at(tt+t);
            factor2 += autocorr_mags.at(tt);
            factor3 += autocorr_mags.at(tt + t);
        }
        ismdl_info.autocorrelation.push_back(f0*factor1 - f0*factor2*f0*factor3);
    }
    
    // normalizes everything so that the autocorrelation at (t=0) = 1
    double a0 = ismdl_info.autocorrelation.at(0);
    for (double &a: ismdl_info.autocorrelation) a /= a0;

    // now we calculate the correlation time
    // we will natural log all of the values,
    // then use GSL to find the best fit. the slope, when inverted and multiplied by -1, is the correlation time
    // NOTE: i have no idea if there is a better way. there is probably some shit buried in <algorithm> that can help, but who cares
    vector<double> autocorrelation_vals;
    for (int i=0; i<CORRELATION_TIME_RANGE; i++) {
        autocorrelation_vals.push_back(log(ismdl_info.autocorrelation.at(i)));
    }

    // need a set of "X" values to do the fit over
    vector<double> N_vals;
    for (int i=0; i<CORRELATION_TIME_RANGE; i++) N_vals.push_back((double)i);

    // compute best fit with GSL
    double b, m, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(N_vals.data(), 1, autocorrelation_vals.data(), 1, CORRELATION_TIME_RANGE, &b, &m, &cov00, &cov01, &cov11, &sumsq);

    // store the slope and y intercept so we can graph them together later
    ismdl_info.auto_correlation_slope = m;
    ismdl_info.auto_correlation_y_int = b;

    cout << m << " ";

    return -1.0/m; // m = -1.0/tau, so tau = -1.0/m.
}   

void IsingModel::update_grid() {

    for (size_t n=0; n<(lattice_size*lattice_size); n++) {
        // pick random lattice point
        int i = rand()%lattice_size, j = rand()%lattice_size;

        // compute potential change in energy
        int delta_U = calc_delta_U(i, j);
        

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


int64_t IsingModel::simulate(int N_iter) {
    auto start = chrono::high_resolution_clock::now();

    for (int i=0; i<N_iter; i++) {
        update_grid();
        ismdl_info.energies.push_back(calc_total_energy_per_site());
        ismdl_info.mags.push_back(calc_total_mag_per_site());
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);

    double correlation_time = calc_autocorrelation();
    cout << "Correlation time: " << correlation_time << "\n";

    return duration.count();
}

int64_t IsingModel::multi_simulate(double T0, double Tf, double dT, double N_iter_equib, double N_iter_avg) {

    auto start = chrono::high_resolution_clock::now();

    int N_iter = (int)((Tf-T0)/dT) + 2; // +1 from double->int, +1 from exclusive upper bound on loop

    for (int n=0; n<N_iter; n++) {

        // calculate and set temperature for this iteration
        double temp = T0 + n*dT;
        set_temp(temp);
        ismdl_info.temps.push_back(temp);

        // first loop gets system in equilibrium
        for (int i=0; i<N_iter_equib; i++) {
            update_grid();
            ismdl_info.intermediate_mags.push_back(abs(calc_total_mag_per_site()));
        }

        double correlation_time = calc_autocorrelation();
        cout << "Correlation time for T=" << temp << ": " << correlation_time << "\n";
        ismdl_info.correlation_time.push_back(correlation_time);

        // second loop is where averages are computed.
        double mag = 0.0;
        double energy = 0.0;
        double mag_squared = 0.0;
        double energy_squared = 0.0;
        for (int i=0; i<N_iter_avg; i++) {
            double _mag = calc_total_mag_per_site();
            double _energy = calc_total_energy_per_site();
            mag += abs(_mag); // absolute value of the mean magnetization for a given run
            energy += _energy;
            mag_squared += _mag*_mag;
            energy_squared += _energy*_energy;

            update_grid();
        }
        // divide the quantities by the number of samples that were taken to get means
        mag /= N_iter_avg;
        energy /= N_iter_avg;
        mag_squared /= N_iter_avg;
        energy_squared /= N_iter_avg;


        // add all of the values to vectors to plot
        ismdl_info.mags.push_back(mag);
        ismdl_info.energies.push_back(energy);
        ismdl_info.specific_heats.push_back(beta*beta*(energy_squared - energy*energy)/lattice_size);
        ismdl_info.susceptibilities.push_back(beta*lattice_size*(mag_squared - mag*mag));


        total_energy = 0; // reset the total energy for the next simulation
        cout << "T=" << temp << " sim complete." << "\n"; // good way to know that it is working!
    }

    auto end = chrono::high_resolution_clock::now();
    cout << "All sims complete\n";

    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start); // prob should return seconds...
    return duration.count();
}


void IsingModel::plots() {
    if (sim_type == SINGLE_SIM) {
        // plotting energy
        plt::title("Energy vs. Iterations per lattice site");
        plt::plot(ismdl_info.energies, "k");
        plt::xlabel("Iterations per lattice site (sweeps) N");
        plt::ylabel("(Dimensionless) Energy E");
        plt::save("res/graphs/single_sims/energy.png");

        plt::cla();

        // plotting magnetization
        plt::title("Magnetization vs. Iterations per lattice site");
        plt::plot(ismdl_info.mags, "k");
        plt::xlabel("Iterations per lattice site (sweeps) N");
        plt::ylabel("Magnetization M");
        plt::save("res/graphs/single_sims/magnetization.png");
        plt::cla();


        // grab only the first few hundred, i.e. the first few correlation times to plot
        // otherwise, long-term exponential behavior would dominate, which we don't want!
        vector<double> autocorrelation_vals(
            ismdl_info.autocorrelation.begin(),
            ismdl_info.autocorrelation.begin() + CORRELATION_TIME_RANGE*4
        );
        // we also need an "X" vector to plot using semilogy. 
        // the docs say I don't need to, but who knows.
        vector<double> N;
        for (size_t i=0; i<CORRELATION_TIME_RANGE*4; i++) N.push_back((double)i);

        plt::title("Autocorrelation vs. Iterations per lattice site");
        plt::plot(N, autocorrelation_vals, "k");
        plt::xlabel("Iterations per lattice site (sweeps) N");
        plt::ylabel("Magnetization autocorrelation $\\chi$");
        plt::save("res/graphs/single_sims/autocorrelation.png");
        plt::cla();



        // we want to choose a smaller set, then plot on logarithmic scale
        vector<double> N_log(N.begin(), N.begin() + CORRELATION_TIME_RANGE*2);
        vector<double> log_autocorrelation(autocorrelation_vals.begin(), autocorrelation_vals.begin() + CORRELATION_TIME_RANGE*2);
        // ln() everything to plot logarithmic
        for_each(log_autocorrelation.begin(), log_autocorrelation.end(), [](double &x){x = log(x); });

        // generate the best fit line
        double m = ismdl_info.auto_correlation_slope;
        double b = ismdl_info.auto_correlation_y_int;
        vector<double> best_fit_line;
        for (double x: N_log) {
            best_fit_line.push_back(m*x + b);
        }

        plt::title("Log Autocorrelation and best fit line");
        plt::plot(N_log, log_autocorrelation, {{"color", "black"},{"label", "Autocorrelation"}});
        plt::plot(N_log, best_fit_line, {{"color", "blue"},{"label", "Best fit line"},{"linestyle", "--"}});
        plt::xlabel("Iterations per lattice site (sweeps) N");
        plt::ylabel("Log Magnetization autocorrelation $\\log \\chi(t)$");
        plt::legend();
        plt::save("res/graphs/single_sims/log_autocorrelation.png");

        plt::detail::_interpreter::kill(); // idk why this is needed on my laptop
    } else if (sim_type == MULTI_SIM) {

        // the padding is to make this look nicer. 
        // if I just limited it to the minimum, half of a dot would be off the screen
        // makes a lot of extra lines and calculations, but compared to the entire multi-simulation,
        // this is nearly nothing.
        double pad = 0.0;

        // magnetization vs. temperature
        double max_mag = *max_element(ismdl_info.mags.begin(), ismdl_info.mags.end());
        double min_mag = *min_element(ismdl_info.mags.begin(), ismdl_info.mags.end());
        double mag_range = max_mag - min_mag;
        pad = mag_range*0.05;
        plt::title("Magnetization vs. Temperature");
        plt::plot(ismdl_info.temps, ismdl_info.mags, "k.");
        plt::plot({T_EXACT,T_EXACT}, {min_mag - pad, max_mag + pad}, {{"linestyle","--"},{"color","blue"}});
        plt::ylim(min_mag - pad, max_mag + pad);
        plt::xlabel("Temperature T");
        plt::ylabel("Magnetization m");
        plt::save("res/graphs/multi_sims/magnetization.png");
        plt::cla();


        // specific heat vs. temperature. this one is fucked, i don't know why!
        // the code to calculate this is identical to what I had a few days ago, 
        // when everything was working fine. who knows!
        double max_specific_heat = *max_element(ismdl_info.specific_heats.begin(), ismdl_info.specific_heats.end());
        double min_specific_heat = *min_element(ismdl_info.specific_heats.begin(), ismdl_info.specific_heats.end());
        double specific_heat_range = max_specific_heat - min_specific_heat;
        pad = specific_heat_range*0.05;
        plt::title("Specific Heat vs. Temperature");
        plt::plot(ismdl_info.temps, ismdl_info.specific_heats, "k.");
        plt::plot({T_EXACT,T_EXACT}, {min_specific_heat-pad,max_specific_heat+pad}, {{"linestyle","--"},{"color","blue"}});
        plt::ylim(min_specific_heat-pad,max_specific_heat+pad);
        plt::xlabel("Temperature T");
        plt::ylabel("Specific Heat c");
        plt::save("res/graphs/multi_sims/specific_heats.png");
        plt::cla();


        // double max_susceptibility = *max_element(ismdl_info.susceptibilities.begin(), ismdl_info.susceptibilities.end());
        double max_susceptibility = 0.15; // hard coded because there is one value that goes super high and messes up scaling
        double min_susceptibility = *min_element(ismdl_info.susceptibilities.begin(), ismdl_info.susceptibilities.end());
        double susceptibility_range = max_susceptibility - min_susceptibility;
        pad = susceptibility_range*0.05;
        plt::title("Susceptibility vs. Temperature");
        plt::plot(ismdl_info.temps, ismdl_info.susceptibilities, "k.");
        plt::plot({T_EXACT,T_EXACT}, {min_susceptibility-pad,max_susceptibility+pad}, {{"linestyle","--"},{"color","blue"}});
        plt::ylim(min_susceptibility-pad,max_susceptibility+pad);
        plt::xlabel("Temperature T");
        plt::ylabel("Susceptibility \\chi");
        plt::save("res/graphs/multi_sims/susceptibility.png");
        plt::cla();

        double max_correlation_time = *max_element(ismdl_info.correlation_time.begin(), ismdl_info.correlation_time.end());
        double min_correlation_time = *min_element(ismdl_info.correlation_time.begin(), ismdl_info.correlation_time.end());
        double correlation_time_range = max_correlation_time - min_correlation_time;
        pad = correlation_time_range*0.05;
        plt::title("Correlation time vs. Temperature");
        plt::plot(ismdl_info.temps, ismdl_info.correlation_time, "k.");
        plt::plot({T_EXACT,T_EXACT}, {min_correlation_time-pad,max_correlation_time+pad}, {{"linestyle","--"},{"color","blue"}});
        plt::ylim(min_correlation_time-pad,max_correlation_time+pad);
        plt::xlabel("Temperature T");
        plt::ylabel("Correlation time \\tau");
        plt::save("res/graphs/multi_sims/correlation_time.png");
        plt::cla();



        plt::detail::_interpreter::kill(); // idk why this is needed on my laptop
    }
}


void IsingModel::draw_grid() {
    int w = raylib_info.win_w;
    int h = raylib_info.win_h;
    int n = lattice_size;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            Color color = (grid[i][j] > 0) ? WHITE : BLACK;
            DrawRectangle(j*(w/n), i*(h/n), (w/n), (h/n), color);
        }
    }
}

void IsingModel::show_simulation() {
    InitWindow(raylib_info.win_w, raylib_info.win_h, raylib_info.win_title.c_str());
    SetTargetFPS(raylib_info.refresh_rate);

    // this is a little silly
    string screenshot_output_path = SCREENSHOT_OUTPUT_PATH;
    screenshot_output_path += "/";

    int iter_counter = 0;
    int screenshot_counter = 0;

    while (!WindowShouldClose()) {
        BeginDrawing();

        // want to be able to restart it if we press "r"
        if (IsKeyPressed(KEY_R)) {
            init_grid(grid_init_flag);
            iter_counter = 0;
            screenshot_counter = 0;
        }

        if (iter_counter == raylib_info.screenshots[screenshot_counter]) {
            string file_name = screenshot_output_path + "screenshot" + to_string(iter_counter) + ".png";
            string file_name2 = "../Report/figures/screenshot" + to_string(iter_counter) + ".png";

            cout << file_name2 << "\n";

            // TakeScreenshot is supposed to be a wrapper around this, but for some reason,
            // it won't accept relative paths...
            // Guess I gotta do it the long way
            Image i = LoadImageFromScreen();
            ExportImage(i, file_name.c_str());
            ExportImage(i, file_name2.c_str()); 
            screenshot_counter++;
        }

        ClearBackground(BLACK);
        draw_grid();
        update_grid();

        iter_counter++;

        EndDrawing();
    }

    CloseWindow();
}




void IsingModel::run() {
    // each part of the if(else) statement just grabs the values specific for that run,
    // then executes the corresponding function

    if (sim_type == SINGLE_SIM) {
        double equib_temp = stod(init_data["EQUIB_TEMP"]);
        set_temp(equib_temp);
    
        // if we have equilibration time override or not
        int equib_override = stoi(init_data["OVERRIDE_EQUIB_TIME"]);
        int N_iter = (equib_override) ? stoi(init_data["EQUIB_TIME"]) : NUM_EQUIB_ITER;

        cout << "Running single-simulation...\n";
        int64_t ms = simulate(N_iter);
        cout << "\nDone! Simulation took " << ms << " milliseconds\n";
        plots();
    } else if (sim_type == MULTI_SIM) {
        double t0 = stod(init_data["T0"]);
        double tf = stod(init_data["TF"]);
        double dt = stod(init_data["DT"]);

        cout << "Running multi-simulation...\n";
        int64_t ms = multi_simulate(t0, tf, dt, NUM_EQUIB_ITER, NUM_AVG_ITER);
        cout << "\nDone! All simulations took " << ms << " milliseconds, or " << (double)ms/1000.0 << " seconds.\n";
        plots();
    } else if (sim_type == RAYLIB_SIM) {

        raylib_info.win_w = stoi(init_data["WIN_W"]);
        raylib_info.win_h = stoi(init_data["WIN_H"]);
        raylib_info.win_title = init_data["WIN_TITLE"];
        raylib_info.refresh_rate = stoi(init_data["REFRESH_RATE"]);

        raylib_info.equib_temp = stod(init_data["EQUIB_TEMP"]);
        set_temp(raylib_info.equib_temp);

        // grab all the screenshots numbers. Even in C++, this is annoying! Could be a one-liner with python...
        string screenshots_str = init_data["SCREENSHOTS"];
        unsigned int comma_pos = 0;
        int counter = 0;
        while (comma_pos < screenshots_str.size()) {
            if (counter == MAX_SCREENSHOTS) break;

            comma_pos = screenshots_str.find(",");
            string s = screenshots_str.substr(0,comma_pos);
            int screenshot = stoi(s);
            raylib_info.screenshots[counter++] = (screenshot == 0) ? 1 : screenshot; //avoid putting zero
            screenshots_str.erase(0,comma_pos+1);
        }

        show_simulation();
    }
}