[![DOI](https://zenodo.org/badge/1052689875.svg)](https://doi.org/10.5281/zenodo.17078682)

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

# pyDRACULA

DRACULA: Decay, Recombination And Collisions compUtation for Low-temperature Antihydrogen

## General:

Python version (for C++ version, see https://github.com/balintradics/Dracula)

This is a code to calculate the formation and level population of antihydrogen atoms, from a single passage of antiprotons through a positron plasma. Starting from empty level population (an ensemble of bare antiprotons) antihydrogen bound states are formed via three-body and radiative recombination processes. The antihydrogen atoms can also be (de)excited, and ionised via collisions with positrons. All these processes are calculated by time evolving the corresponding coupled rate equations for antiproton and antihydrogen levels. For further details please see the paper below.

If you use the code please cite: Phys. Rev. A 90, 032704 (2014), B. Radics, D.J. Murtagh, Y. Yamazaki and F. Robicheaux
DOI:https://doi.org/10.1103/PhysRevA.90.032704

## Requirements:

sqlite3 (tested 3.20.1), python (tested with 2.7), pygsl (2.2.0), matplotlib, numpy

## Usage example:

To use the code, first unzip the file `hydrogenDB.zip` and name it to `hydrogenDB.db`.
The code needs to be steered by a python script. Here is a self-explanatory example:

```
import matplotlib.pyplot as plt
from numpy import array
import numpy as np
from Dracula import Dracula


# 1. Overlap a dense positron plasma and a bare antiproton plasma
# Units: positron_temp (K), positron_density (m^{-3}), antiproton_number (1), magnetic_field (T), time (microseconds), radiationTemp (K)

d = Dracula(NStates = 150, positron_temp =50, positron_density = 1e14, antiproton_number = 100000, time = 20, magnetic_field= 2, radiationTemp = 300, useBlackBody = True, saveGSSteps = True)

# Read the various rate coefficients
d.SetSQLDbScat("hydrogenDB.db", True)
d.SetSQLDbRadDec("hydrogenDB.db", True)

# Evolve the system to form bound states
d.Compute()

# Save the level population into CSV file (Principal QN, Population)
# Output level population unit is absolute number of antihydrogens
d.SaveLevelPopCSV("overlap.csv")

# Save the states and steps for further processing
levpop = d.GetLevelPop()
steps = d.GetSaveSteps()
steps = array(steps)
```

The running and the output of this code (see also in `main1.py`) should be something similar as below:

```
$ python3.9 main1.py
Reading in scattering rates from sqlite3 file hydrogenDB.db ...
Reading radiative recombination rates from sqlite3 file hydrogenDB.db...
Reading in decay rates from sqlite3 file hydrogenDB.db...this might take a while...
Beginning computation...
Initial level population...
[     0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0.      0.      0.      0.
      0.      0.      0.      0.      0.      0. 100000.]
# Using stepper bsimp with order 12
Starting integration...
Time, Stepsize, GS lev.pop.:  1e-09 5e-09 2.6180173762271302e-08
Time, Stepsize, GS lev.pop.:  6e-09 2.0098642105860833e-08 2.071596291781851e-07
Time, Stepsize, GS lev.pop.:  2.6098642105860834e-08 6.445684635190868e-08 1.2080067941320795e-06
Time, Stepsize, GS lev.pop.:  9.055548845776951e-08 1.9283124624884672e-07 6.177346616506293e-06
Time, Stepsize, GS lev.pop.:  2.8338673470661623e-07 5.367713321628146e-07 3.920845870959867e-05
Time, Stepsize, GS lev.pop.:  8.201580668694309e-07 1.432841847333307e-06 0.0002646255156388598
Time, Stepsize, GS lev.pop.:  2.252999914202738e-06 3.5232744129159972e-06 0.001439768926144104
Time, Stepsize, GS lev.pop.:  5.776274327118735e-06 8.979483285975207e-06 0.006943831374921808
Time, Stepsize, GS lev.pop.:  1.4755757613093942e-05 2.3365864256816616e-05 0.043651905234487906
Time, Stepsize, GS lev.pop.:  1.9999999999999998e-05 2.622121193453028e-05 0.08559158162159901
Total number of bound states :20.075127807168034, number of GS antiH: 0.08559158162159901, burned antiprotons 20.075127806601813
Saving result to csv file overlap.csv
```
The contents (principal QN, population) of the output CSV can be easily plotted:

```
import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('overlap.csv', delimiter=',', skip_header=0,
                     skip_footer=0, names=['n', 'h'])

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.semilogy(data['n'], data['h'], color='r', marker='.', linestyle='None', label='overlap')
ax1.set_xlabel('N principal qn')
ax1.set_ylabel('N. or atoms [a.u.]')

plt.legend(bbox_to_anchor=(0.8, 0.5), loc=2, borderaxespad=0.)

plt.grid(True)

plt.show()
```
This should result in a plot similar to the following:

![Antihydrogen level population](/Images/Figure_1.png)
