import matplotlib.pyplot as plt
from numpy import array
import numpy as np
from Dracula import Dracula


# 1. Overlap a dense positron plasma and a bare antiproton plasma 

d = Dracula(NStates = 150, positron_temp =50, positron_density = 1e14, antiproton_number = 100000, time = 20, magnetic_field= 2, radiationTemp = 300, useBlackBody = True, saveGSSteps = True)

# Read the various rate coefficients
d.SetSQLDbScat("hydrogenDB.db", True)
d.SetSQLDbRadDec("hydrogenDB.db", True)

# Evolve the system to form bound states
d.Compute()

# Save the level population into CSV file (Principal QN, Population)
d.SaveLevelPopCSV("overlap.csv")

# Save the states and steps for further processing
levpop = d.GetLevelPop()
steps = d.GetSaveSteps()
steps = array(steps)


# 2. Use the previous level population output,
# but allow only in-flight rad. decay, turn off all scattering processes
# so antiproton number equals to 0

d2 = Dracula(NStates = 150, positron_temp =50, positron_density = 1e14, antiproton_number = 0, time = 100, magnetic_field= 2, radiationTemp = 300, useBlackBody = True, useInitialLevelPop = True, saveGSSteps = True)

# Read the various rate coefficients
d2.SetPop(levpop)
d2.SetSQLDbRadDec("hydrogenDB.db", True)

# Evolve the system to form bound states
d2.SetStartTime(steps[-1,:-1])
d2.Compute()

# Save the level population into CSV file (Principal QN, Population)
d2.SaveLevelPopCSV("flight.csv")

# Save the states and steps for further processing
levpop2 = d2.GetLevelPop()
steps2 = d2.GetSaveSteps()
steps2 = array(steps2)

stepsall = np.concatenate((steps, steps2))

print("Time steps [s], GS lev.pop.")
print(stepsall)

