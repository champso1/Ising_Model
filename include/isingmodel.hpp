#pragma once

#include <string_view>
#include <vector>
#include <random>
#include <iostream>

namespace ising
{
	// Default simulation values
	inline constexpr size_t num_equib_iter = 5000;
	inline constexpr size_t num_avg_iter = 15000;

	// Other default vals for calculations
	inline constexpr double t_exact = 2.269;

	inline std::mt19937 gen(std::random_device{}());
	// returns a random number between 0 and 1
	double rand_double();



	inline constexpr char const* ANSI_COLOR_RED =         "\x1b[31m";
	inline constexpr char const* ANSI_COLOR_GREEN =       "\x1b[32m";
	inline constexpr char const* ANSI_COLOR_YELLOW =      "\x1b[33m";
	inline constexpr char const* ANSI_COLOR_BLUE =        "\x1b[34m";
	inline constexpr char const* ANSI_COLOR_MAGENTA =     "\x1b[35m";
	inline constexpr char const* ANSI_COLOR_CYAN =        "\x1b[36m";
	inline constexpr char const* ANSI_COLOR_RESET =       "\x1b[0m";
	inline constexpr char const* ANSI_LINEFEED_CLEAR =    "\033[2K";
	inline constexpr std::string_view loading_block("█");

	void log_iterations(size_t val, size_t end, std::string_view name);
	inline void finish_log_iterations()
	{
		std::cout << std::endl << std::flush;
	}
	


	// Grid init flags
	enum class GridInitFlag
	{
		INIT_T_0, // all spins one direction
		INIT_T_INF // all spins completely random
	};

	// quick helper class to keep grid data in one flat vector
	// supplying a getter to handle strides correctly
	class Grid final
	{
	public:
		using grid_type = std::vector<int>;
	private:
		size_t _rows, _cols;
		grid_type _grid;
	public:
		Grid() = delete;
		Grid(size_t rows, size_t cols);
		void fill(GridInitFlag init_flag);
		inline int const& at(size_t row, size_t col) const { return _grid[row*_cols+col]; }
		inline int& at(size_t row, size_t col) { return _grid[row*_cols+col]; }
		inline auto size() const { return _rows; }
	};


	class IsingModel final
	{
	public:
		using sim_return_type = std::pair<std::vector<double>,std::vector<double>>;
		
		struct MultisimResult final
		{
			double avg_energy;
			double avg_magnetization;
			double susceptibility;
			double specific_heat;
		};
		using multisim_return_type = std::vector<std::pair<double,MultisimResult>>;
	private:
		Grid _grid;
	    GridInitFlag _grid_init_flag;

		double _beta; // 1/T
		int _total_energy; // current total energy in the lattice
	public:
		IsingModel() = delete;
		// initializes the grid and everything
		IsingModel(size_t N, GridInitFlag grid_init_flag);
		~IsingModel() = default;
		
		// runs one simulation at one temperature
		// returns number of seconds it took
		sim_return_type simulate(double equil_temp, int N_iter);

		// runs multiple simulations in a range of temperatures
		multisim_return_type multi_simulate(
			double T0,
			double Tf,
			double dT,
			double N_iter_equib,
			double N_iter_avg);

#if ENABLE_RAYLIB
		// uses raylib to show a single simulation run
		// doesn't do any calculations
		void show_simulation(
			double equil_temp,
			int win_w, int win_h,
			char const* win_title,
			int refresh_rate
		);
#endif

	private:
		// setups up the grid either all zeros/ones or all random
		void init_grid(GridInitFlag flag);

		// calculates the change in energy due to flipping one spin
		// used to determine if a proposed metropolis step should be accepted
		int calc_deltaE(int i, int j); 
		// performs a full metropolis run, updating the entire grid
		void update_grid();

		// calculates the total energy of the lattice
		double calc_energy();
		// calculates the total magnetization of the whole lattice
		double calc_magnetization();

		// some helpers
		inline void set_temp(double temp) { _beta = 1.0/temp; }
		inline void reset_energy() { _total_energy = 0; } // for clarity inside of functions
	};
}
