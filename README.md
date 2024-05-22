# Two-Dimensional Ising Model Simulation
This project was done for my PHYS3500K Computational Physics I project.

## Computational Methods
This simulation relies on the Algorithm of Metropolis et. al. The report explaining the details can be found in the Report which I have added to this repository, but a quick recap is as follows:

The Metropolis algorithm specializes in generating random values according to specified probability distributions, and it is often very good at doing so; far better than other naive methods, despite it being pretty simple. It achieves this by sampling a particlar value in the domain of the problem, then generating a proposal sample based on the previous one. It will add this sample to its distribution based on a weight function characteristic of the probability distribution you are attempting to sample. 

In the case of the Ising model, this weight function is the ratio of the Boltzmann factors of the two states, which we are allowed to use because we are working in the canonical ensemble (we provide the temperature, the number of dipoles, and the area) It can be shown that by doing this, the algorithm produces subsequent states that satisfy the Boltzmann distribution, which is exactly what we want.

From here, we can start an NxN lattice off at some initial temperature, say infinity or zero (corresponding to dipoles being randomly spin up/down, or all in one direction, respectively), then let it come to equilibrium at some other specified temperature. We can perform multiple of these simulations in a row, letting it come to equilibrium at various temperatures in order to probe quantities such as the magnetization, specific heat, susceptibility, or the autocorrelation, and examine what happens at the critical point, found analytically by Lars Onsager in 1944.


## Installation and Use
This project requires a few dependencies, and it is highly recommended that this is run on a Unix-based system (MacOX or Linux). Ensure that you have the pkg-config package installed on your system BEFORE doing this, as the Makefile relies on them, and if it is not installed, neither Raylib nor GSL will install a `.pc` file. This can be done with `sudo apt-get install pkg-config` on Ubuntu or `brew install pkg-config` on MacOS, for instance.

The main dependencies are [Raylib](https://github.com/raysan5/raylib) for the graphical simulation, [GSL](https://www.gnu.org/software/gsl/) for an efficient least-squares fitting algorithm which is used to find the correlation time, and [ROOT](https://root.cern/),for data storage and plotting. Their corresponding pages have detailed installation instructions.


Once these are all installed, run the following to do the simulation:

```
git clone https://github.com/champso1/Ising_Model.git
cd Ising_Model
make
./bin/isingmodel
```

after setting details inside of the `run.config` file.

To plot the results after a simulation has been run and a `.root` file has been created with the resultant data, you can run:

```
root ./res/plot.C
```

from the terminal, or 

```
.x ./res/plot.C
```

from inside the ROOT shell.