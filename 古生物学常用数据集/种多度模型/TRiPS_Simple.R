# Script to perform estimation of sampling rates and species richness using TRiPS.
# Written by 
# Jostein Starrfelt (jostein.starrfelt@ibv.uio.no), 
# Centre for Ecological and Evolutionary Synthesis, University of Oslo, Norway.

# The code is an example on how to apply TRiPS to a set of occurrence counts inside
# a geological interval of a given duration (dt). Please see main manuscript for details
# and caveats.
# Main manuscript.
# "How many dinosaur species were there? True richness estimated using a Poisson sampling model (TRiPS)"
# Authors: Jostein Starrfelt & Lee Hsiang Liow.
# To be printed in Philosophical Transactions B. 

# stats4 library is required. If you do now have it, install the package before
# running this code.
library(stats4)
# The file below has all the necesary function for TRiPS.
source("functions_DinoBino_toweb.R")


# TRiPS needs two kinds of input; a set of occurrence counts (a list of how many occurrences
# per species has been observed) and the duration of the interval in which they were
# observed.

# Say we have observed 5 species with occurrence counts (4,3,1,1,1) in an interval of duration
# 10 Ma.
obs = c(4,3,1,1,1)  # occurrence counts
n_raw = length(obs) # Observed richness
dt  = 10            # Duration of interval.

# All calculations have been subsumed in the function doTRiPS_abs which takes
# the occurrence count lists and the duration of the interval as inputs and outputs
# maximum likelihood estimates for sampling rate, sampling probability and the estimated
# true richness.

Out = doTRiPS_abs(obs,t=10)

# End