# Two-Dimensional Ising Model Simulation
This project was done for my PHYS3500K Computational Physics I project.

## Computational Methods
This simulation relies on the Algorithm of Metropolis et. al. The report explaining the details can be found in the Report which I have added to this repository, but a quick recap is as follows:

The Metropolis algorithm specializes in generating random values according to specified probability distributions, and it is often very good at doing so; far better than other naive methods, despite it being pretty simple. It achieves this by sampling a particlar value in the domain of the problem, then generating a proposal sample based on the previous one. It will add this sample to its distribution based on a weight function characteristic of the probability distribution you are attempting to sample. 

In the case of the Ising model, this weight function is the ratio of the Boltzmann factors of the two states, which we are allowed to use because we are working in the canonical ensemble (we provide the temperature, the number of dipoles, and the area) It can be shown that by doing this, the algorithm produces subsequent states that satisfy the Boltzmann distribution, which is exactly what we want.

From here, we can start an NxN lattice off at some initial temperature, say infinity or zero (corresponding to dipoles being randomly spin up/down, or all in one direction, respectively), then let it come to equilibrium at some other specified temperature. We can perform multiple of these simulations in a row, letting it come to equilibrium at various temperatures in order to probe quantities such as the magnetization, specific heat, susceptibility, or the autocorrelation, and examine what happens at the critical point, found analytically by Lars Onsager in 1944.


## Installation and Use
This project requires a few dependencies, and it is highly recommended that this is run on a Unix-based system (MacOX or Linux). Ensure that you have the pkg-config package installed on your system BEFORE doing this, as the Makefile relies on them, and if it is not installed, neither Raylib nor GSL will install a `.pc` file. This can be done with `sudo apt-get install pkg-config` on Ubuntu or `brew install pkg-config` on MacOS, for instance.

The main dependencies are [Raylib](https://github.com/raysan5/raylib) for the graphical simulation and [GSL](https://www.gnu.org/software/gsl/) for an efficient least-squares fitting algorithm which is used to find the correlation time. Their corresponding pages have detailed installation instructions. Raylib requires a few other graphic API dependences, but can be build easily with CMake. If you build shared libraries for Raylib (as I have done), ensure you add `/usr/local/lib` to your `LD_LIBRARY_PATH` if you are working on Linux. This is a common problem, and can be found at the end of the corresponding [wiki page](https://github.com/raysan5/raylib/wiki/Working-on-GNU-Linux). GSL has a `./configure; make; make install;` system that is very simple.

Additionally, ensure you have Python development libraries. For instance, on Ubuntu 22.04, ensure you `sudo apt-get install python3-dev python3-numpy python3-matplotlib` to get development headers for those libraries. I use [Matplotlibcpp](https://github.com/lava/matplotlib-cpp) to pipe plotting data to a Python interpreter equipped with Matplotlib, as that is the syntax I am familiar with, and it requires these headers.

Once these are all installed, run the following:

```
git clone LINK
cd LINK
make
./bin/isingmodel
```

After you have set up the `run.config` file to your liking.