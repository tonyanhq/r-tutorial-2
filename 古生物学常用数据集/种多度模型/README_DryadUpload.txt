README for files containing R-code and data for TRiPS manuscript;
"How many dinosaur species were there? Fossil bias and true richness estimated using a
Poisson sampling model." 
by Jostein Starrfelt & Lee Hsiang Liow
email: jostein.starrfelt@ibv.uio.no
Phil Trans B 2016, doi: 10.1098/rstb.2015.0219.

This upload contains several scripts, functions and RData files with results and algorithms for TRiPS analysis on dinosaur data as well as simulations used in the main manuscript. Code is admittedly not optimized, but supplied for interested parties to replicate the analyses. Code was generated in Rstudio, but should work for most R applications.

This code is licensed under the CC0 license. Do whatever you wish, as long as the original author is properly attributed.

-- Simple application of TRiPS --

To apply TRiPS to estimate sampling rates and probabilities as well as true richness the R script TRiPS_Simple.R shows a simple example.


-- Main analysis in Phil Trans --

This part has 5 .R files where TRiPS_MainsScript_Revision.R details all aspects of analysis, while functions_DinoBino_toweb.R contains the necessary functions and 3 scripts for plotting the figures.

To reproduce whole or part of the analysis the first part of TRiPS_MainScript has a series of switches. Setting them as TRUE/FALSE will in various degrees replicate our analysis when running the whole script. The switches are set to load the results and plot the figures, but can also be set to re-run analysis on current data, or to download new data. This requires a few other R packages, which are obvious from the initial few lines.

Scripts/functions:
TRiPS_MainScript_Revision.R
functions_DinoBino_toweb.R
Figure_1_Sampling_191115.R
RichnessFig_rev.R
plotallsamplingrates_SI.R

- Data:
Raw: 
PBDB_Download_only.RData
invalid_PBDB.RData
Ootaxa_Benson.txt
Ichnotaxa_Benson.txt

Results:
countedcollections.RData
RevisedAnalysis_Genera.RData
RevisedAnalysis_Species.RData

Please note that there are also a few Oo- and ichno-taxa manually removed inside the script, as they were not in the lists supplied by Roger Benson (see main script).

-- Simulations --

This part has one script and one RData file with stored simulation results. The script has one switch (dosimulation) which as default has FALSE, which when set to TRUE will make the script re-run simulations, and then plot new figures. Number of simulations should then be set (in the manuscript we used 100000). Parameter ranges could also be changed [lines 22-44 in SimulationScript_withFigures.R]

Scripts/functions:
SimulationScript_withfigures.R

Data/output from simulations used in manuscript:
Simus_run_220915.RData




All files uploaded:
Figure_1_Sampling_191115.R
functions_DinoBino_toweb.R
plotallsamplingrates_SI.R
RichnessFig_rev.R
Simulate_BDF_functRat_newLvar.R
SimulationScript_withfigures.R
TRiPS_MainScript_Revision.R
PBDB_Download_only.RData
invalid_PBDB.RData
countedcollections.Rdata
RevisedAnalysis_Genera.RData
RevisedAnalysis_Species.RData
Simus_run_220915.RData
README_DryadUpload.txt
Ichnotaxa_Benson.txt
Ootaxa_Benson.txt
