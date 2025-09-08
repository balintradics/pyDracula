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



