#include "isingmodel.hpp"
using namespace ising;

#include <fstream>

static void output_singlesim_results(IsingModel::sim_return_type const& singlesim_res);

int main()
{
	IsingModel ismdl(400, GridInitFlag::INIT_T_INF);
	
	IsingModel::sim_return_type res = ismdl.simulate(3.0, 1000);
	output_singlesim_results(res);
	
    return 0;
}


static void output_singlesim_results(IsingModel::sim_return_type const& singlesim_res)
{
	auto [energies,mags] = singlesim_res;
	std::ofstream outfile("single.dat");
	for (size_t i=0; i<energies.size(); ++i)
		outfile << i << ' ' << energies[i] << ' ' << mags[i] << '\n';
		
}
