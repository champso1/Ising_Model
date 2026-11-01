#include "isingmodel.hpp"

#include <algorithm>
#include <print>
#include <chrono>
#include <ranges>
#include <cmath>

namespace ising
{
	void log_iterations(size_t val, size_t end, std::string_view name)
	{
	    double ratio = static_cast<double>(val)/static_cast<double>(end);
		int num_blocks = static_cast<int>(ratio*50.0);
		std::string blocks{};
		for (int i=0; i<num_blocks; ++i)
			blocks += loading_block;

		std::string all_text = std::format("\r{}: Iteration {:0>2}/{} ({: >3}%) [{: <50}]",
			name,
			val, end, static_cast<size_t>(ratio*100.0),
			blocks);
		std::cout << all_text << std::flush;
	}


	
	Grid::Grid(size_t rows, size_t cols)
		: _rows{rows}, _cols{cols},
		  _grid(rows*cols, 0.0)
	{}

	void Grid::fill(GridInitFlag init_flag)
	{
		if (init_flag == GridInitFlag::INIT_T_0) {
			std::ranges::fill(_grid, 1);
		} else if (init_flag == GridInitFlag::INIT_T_INF) {
			std::ranges::generate(_grid, [](){ return rand_double() > 0.5 ? 1 : -1; });
		}
	}
	
	IsingModel::IsingModel(size_t N, GridInitFlag grid_init_flag)
		: _grid(N, N), _grid_init_flag{grid_init_flag}
	{
		_grid.fill(grid_init_flag);
	}


	double IsingModel::calc_energy()
	{
		auto lattice_size = _grid.size();
		
		if (_total_energy != 0)
			return static_cast<double>(_total_energy)/(lattice_size*lattice_size); // it will have already been calculated

		// otherwise, we need to calculate it! this only is reached in a call to init_grid()
		int energy = 0;
		for (size_t i=0; i<lattice_size; i++) {
			for (size_t j=0; j<lattice_size; j++) {
				int i_minus1 = (i+lattice_size-1)%lattice_size;
				int j_minus1 = (j+lattice_size-1)%lattice_size;
				int i_plus1 = (i+1)%lattice_size;
				int j_plus1 = (j+1)%lattice_size;
            
				int energy_neighbors = -(_grid.at(i_plus1,j) + _grid.at(i_minus1,j) + _grid.at(i,j_plus1) + _grid.at(i,j_minus1));
				int energy_current = _grid.at(i,j)*energy_neighbors;

				energy += energy_current;
			}
		}
		_total_energy = energy/2; // the above method double counts, so we fix that by dividing by 2

		return static_cast<double>(_total_energy)/(lattice_size*lattice_size);
	}

	double IsingModel::calc_magnetization()
	{
		int mag = 0;
		for (size_t i=0; i<_grid.size(); i++) {
			for (size_t j=0; j<_grid.size(); j++) {
				mag += _grid.at(i,j);
			}
		}
		return (double)mag/(_grid.size()*_grid.size());
	}

	int IsingModel::calc_deltaE(int i, int j)
	{
		auto lattice_size = _grid.size();
		
		int i_minus1 = (i+lattice_size-1)%lattice_size;
		int j_minus1 = (j+lattice_size-1)%lattice_size;
		int i_plus1 = (i+1)%lattice_size;
		int j_plus1 = (j+1)%lattice_size;
    
		int energy_neighbors = -(_grid.at(i_plus1,j) + _grid.at(i_minus1,j) + _grid.at(i,j_plus1) + _grid.at(i,j_minus1));
		int energy_current = _grid.at(i,j)*energy_neighbors;
		int energy_proposed = -energy_current;
		int delta_U = energy_proposed - energy_current;

		return delta_U;
	}

	void IsingModel::update_grid()
	{
		auto lattice_size = _grid.size();
		auto get_rand_index = [lattice_size]() {
			return std::clamp(
				size_t{0}, lattice_size,
				static_cast<size_t>(std::floor(lattice_size*rand_double())));
		};
		
		for (size_t n=0; n<(lattice_size*lattice_size); n++) {
			// pick random lattice point
			
			int i = get_rand_index();
			int j = get_rand_index();

			// compute potential change in energy
			int delta_U = calc_deltaE(i, j);

			if (delta_U <= 0) {         // if change results in lower (or zero) energy, flip automatically
				_grid.at(i,j) *= -1;
				_total_energy += delta_U;
			} else {                    // otherwise, we use the exponential Boltzmann factor as a probability to accept anyway
				double fact = exp(-_beta*delta_U);
				double r = rand_double();
				if (r < fact) {
					_grid.at(i,j) *= -1;
					_total_energy += delta_U;
				}
			}
		}
	}

	IsingModel::sim_return_type IsingModel::simulate(double equil_temp, int N_iter) {
		set_temp(equil_temp);
		
		auto t0 = std::chrono::high_resolution_clock::now();

		std::vector<double> energies(N_iter, 0.0);
		std::vector<double> mags(N_iter, 0.0);
		for (int i=0; i<N_iter; i++) {
			log_iterations(i+1, N_iter, "single sim");
			update_grid();
			energies[i] = calc_energy();
			mags[i] = calc_magnetization();
		}
		finish_log_iterations();

		auto tf = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::ratio<1>> secs = tf-t0;

		std::println("simulation took: {} seconds", secs.count());

		return std::make_pair(energies, mags);
	}
	
	IsingModel::multisim_return_type IsingModel::multi_simulate(
		double T0, double Tf, double dT,
		double N_iter_equib, double N_iter_avg)
	{
		auto t0 = std::chrono::high_resolution_clock::now();

		namespace rv = std::views;
		size_t num_temps = (Tf-T0)/dT + 1;
		auto temps =
			rv::iota(size_t{0}, num_temps)
			| rv::transform([&](size_t i){ return T0 + static_cast<double>(i)*(Tf-T0)/(num_temps-1); });

		multisim_return_type res{};
		for (auto temp : temps) {

			// calculate and set temperature for this iteration
			set_temp(temp);

			// first loop gets system in equilibrium
			for ([[maybe_unused]] int i=0; i<N_iter_equib; i++)
				update_grid();

			// second loop is where values/averages are computed.
			double avg_mag = 0.0;
			double avg_energy = 0.0;
			double mag_squared = 0.0;
			double energy_squared = 0.0;
			for (int i=0; i<N_iter_avg; i++) {
				double energy = calc_energy();
				double mag = calc_magnetization();
				avg_mag += std::abs(mag); // absolute value of the mean magnetization for a given run
				avg_energy += energy;
				mag_squared += mag*mag;
				energy_squared += energy*energy;

				update_grid();
			}
			// divide the quantities by the number of samples that were taken to get means
			avg_mag /= N_iter_avg;
			avg_energy /= N_iter_avg;
			mag_squared /= N_iter_avg;
			energy_squared /= N_iter_avg;

			double suscept = _beta*_grid.size()*(mag_squared - avg_mag*avg_mag);
			double heat = _beta*_beta*(energy_squared - avg_energy*avg_energy)/_grid.size();

			res.emplace_back(std::make_pair(temp, MultisimResult{avg_energy, avg_mag, suscept, heat}));

			reset_energy(); // reset the total energy for the next simulation
			std::println("T = {} simulation complete.", temp);
		}

		auto tf = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::ratio<1>> secs = tf-t0;
		std::println("All sims complete, took {} seconds", secs.count());

		return res;
	}
}
