#include "isingmodel.hpp"

#include <random>

namespace ising
{
	double rand_double()
	{
		static std::uniform_real_distribution<> dist{};
		return dist(gen);
	}

}
