# Script to perform estimation of sampling rates and past species richness of dinosaurs.
# Written by 
# Jostein Starrfelt (jostein.starrfelt@ibv.uio.no), 
# Centre for Ecological and Evolutionary Synthesis, University of Oslo, Norway.

# The code below was used for the analysis presented in
# "How many dinosaur species were there? Fossil bias and true richness estimated using a Poisson sampling model (TRiPS)"
# Authors: Jostein Starrfelt & Lee Hsiang Liow.
# To be printed in Philosophical Transactions B. 

# If you are interested in running TRiPS on own data, a separate script details how to
# best proceed, see TRiPS_Simple.R

# This code is licensed under the CC0 license.

## -- Generating Occurrence count matrices and species count matrices. This is done norep times and 
# the mode of the occ.count sets are used, due to random placement of occurrence IF spanning more than
# one interval. 


## =========   Creating main arrays, downloading data and cleaning data =====
## 


rm(list=ls()) # Empty workspace.
set.seed(sqrt(190480)) # Set seed for reproducability. 

# Switches for reproducing analysis in Starrfelt & Liow.

downloadnew = FALSE # if true all dinosaur data are downloaded again from PBDB. if false data are loaded from file. The loaded data is uncleaned and unstructured.
docleaning  = FALSE # If this is TRUE cleaning of the download will be performed. If doanalysis =TRUE this must also be true to clean the raw download from PBDB
doanalysis  = FALSE # if true analysis is rerun. If False output from saved analysis is loaded and plotted
dogenera    = FALSE # is TRUE genus level analysis/plotting is performed. If FALSE analysis is done on species level.
dofigures   = TRUE # if true figures will be generated.

## -- Loading needed libraries.
library(paleobioDB)
library(stats4)
library(RColorBrewer)
source("functions_DinoBino_toweb.R")

## -- Libraries loaded

# Setting number of replicate runs.
norep = 100 # Number of replications for each occurrence count matrix.

if (downloadnew==FALSE){
  # If not doing new PBDB download
  load("PBDB_Download_only.RData") # Raw data downloads (15th August 2015)
  load("invalid_PBDB.RData")       # Invalid taxa (downloaded Nov 2015)
  load("countedcollections.RData")
} else {
  # Forcing cleaning of data
  docleaning = TRUE
  
  ## -- Loading in data from PBDB. This takes a while.
  # This uses the paleobioDB package. Data has been save to disk and loaded above.
  dinos  <- pbdb_occurrences(limit="all", base_name="Dinosauria", show=c("phylo", "ident", "time"))
  ornits <- pbdb_occurrences(limit="all", base_name="Ornithischia", show=c("phylo", "ident", "time"))
  sauros <- pbdb_occurrences(limit="all", base_name="Sauropodomorpha", show=c("phylo", "ident", "time"))
  theros <- pbdb_occurrences(limit="all", base_name="Theropoda", show=c("phylo", "ident", "time"))
  # There is a bug in the paleobioDB package here, it seems like this doubles all occurrences. 
  # At least the occurrance count is twice the one reported in the Navigator online. 
  length(unique(dinos$oid)) # the number of unique occurrence ID's is half the total downloaded.
  # Bug reported to PBDB and maintainer of R-package.
  # Code to remove dulicates.
  dinos_old = dinos;
  ornits_old = ornits;
  sauros_old = sauros;
  theros_old = theros;
  tmpdins = which(duplicated(dinos_old$oid))
  dinos = dinos_old[-tmpdins,]
  tmpdins = which(duplicated(sauros_old$oid))
  sauros = sauros_old[-tmpdins,]
  tmpdins = which(duplicated(ornits_old$oid))
  ornits = ornits_old[-tmpdins,]
  tmpdins = which(duplicated(theros_old$oid))
  theros = theros_old[-tmpdins,]
  
  # Loading list of invalid taxa.
  invaliddinos <- pbdb_taxa(limit=2000,base_name="Dinosauria",vocab="pbdb",status="invalid",show=c("attr","size","phylo"))
  
  # Collection counts
  tmp = dinos$ein>111 & dinos$lin<139; # collection inside the stages we want.
  collections = unique(dinos[tmp,]$cid)
  nocol = array(0,c(27,1))
  for (ii in 1:length(collections)){
    tmp = which(dinos$cid==collections[ii])
    nocol[(max(112,dinos[tmp[1],]$lin)-111):(max(112,dinos[tmp[1],]$ein)-111)]=nocol[(max(112,dinos[tmp[1],]$lin)-111):(max(112,dinos[tmp[1],]$ein)-111)]+1
  } 
  ## -- done loading in data and removing duplicates
}


if (docleaning==TRUE){
  ## Cleaning up the data. 
  # Downloading list of invalid dinosaur taxa
  invalidnames <- invaliddinos[invaliddinos$rank=="species",5]; # Invalid species names. 
  invalidstatus <- invaliddinos[invaliddinos$rank=="species",7]; # What kind of invalid name.
  
  # === Cleaning up data ===
  # This part is repeated 4 times, for Dinosaurs first and then the subclades.
  #
  # == Dinosauria ==
  # Removing occurrences not identified to species level (mra==3)
  which_not_species=which(dinos$mra!=3)
  dinos = dinos[-which_not_species,]; # Removing all occurrences that are not to matched rank of species
  dim(dinos)
  removedoccs=0;
  invalidwhy = NA; # Storing the index into invalidnames for status of the removed invalid taxonomic names.
  iwtx =1;         # Counter for storing the status of invalid taxonomic names.
  
  # Removing invalid/dubious taxa
  for (ii in 1:length(invalidnames)){
    # For each member of the list of invalid dinosaur taxa
    remove = which(grepl(invalidnames[ii],dinos$mna))
    if (length(remove)>0){
      dinos = dinos[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
      invalidwhy[iwtx] = ii
      iwtx = iwtx+1;
    }
  }
  
  # Ootaxa and Ichnotaxa courtesy of Roger Benson (PLoS in rev). This paper is now published, and the 
  # published taxa lists from the Benson paper was updated with our identification of the additional 4 
  # ichnotaxa.
  ootaxa <- read.table('Ootaxa_Benson.txt')
  ichnotaxa <- read.table('Ichnotaxa_Benson.txt')
  # These are lists of genus, families and orders. We will compare all with genus names in the
  # downloaded data.
  
  for (ii in 1:nrow(ootaxa)){
    remove = which(grepl(ootaxa[ii,],dinos$gnl))
    if (length(remove)>0){
      dinos = dinos[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  dim(dinos) # From 6627 to 6288
  
  for (ii in 1:nrow(ichnotaxa)){
    remove = which(grepl(ichnotaxa[ii,],dinos$gnl))
    if (length(remove)>0){
      dinos = dinos[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  dim(dinos)
  
  # == Ornithischians ==
  which_not_species=which(ornits$mra!=3)
  ornits = ornits[-which_not_species,]; # Removing all occurrences that are not to matched rank of species
  dim(ornits)
  removedoccs=0;
  # Removing invalid/dubious taxa
  for (ii in 1:length(invalidnames)){
    # For each member of the list of invalid dinosaur taxa
    remove = which(grepl(invalidnames[ii],ornits$mna))
    if (length(remove)>0){
      ornits = ornits[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  
  for (ii in 1:nrow(ootaxa)){
    remove = which(grepl(ootaxa[ii,],ornits$gnl))
    if (length(remove)>0){
      ornits = ornits[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  for (ii in 1:nrow(ichnotaxa)){
    remove = which(grepl(ichnotaxa[ii,],ornits$gnl))
    if (length(remove)>0){
      ornits = ornits[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  # From 1535 occurrences of species to 1174
  
  # == Sauropodomorpha ==
  which_not_species=which(sauros$mra!=3)
  sauros = sauros[-which_not_species,]; # Removing all occurrences that are not to matched rank of species
  dim(sauros)
  removedoccs=0;
  # Removing invalid/dubious taxa
  for (ii in 1:length(invalidnames)){
    # For each member of the list of invalid dinosaur taxa
    remove = which(grepl(invalidnames[ii],sauros$mna))
    if (length(remove)>0){
      sauros = sauros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  
  for (ii in 1:nrow(ootaxa)){
    remove = which(grepl(ootaxa[ii,],sauros$gnl))
    if (length(remove)>0){
      sauros = sauros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  for (ii in 1:nrow(ichnotaxa)){
    remove = which(grepl(ichnotaxa[ii,],sauros$gnl))
    if (length(remove)>0){
      sauros = sauros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  
  # == Theropoda ==
  which_not_species=which(theros$mra!=3)
  theros = theros[-which_not_species,]; # Removing all occurrences that are not to matched rank of species
  dim(theros)
  removedoccs=0;
  # Removing invalid/dubious taxa
  for (ii in 1:length(invalidnames)){
    # For each member of the list of invalid dinosaur taxa
    remove = which(grepl(invalidnames[ii],theros$mna))
    if (length(remove)>0){
      theros = theros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  
  for (ii in 1:nrow(ootaxa)){
    remove = which(grepl(ootaxa[ii,],theros$gnl))
    if (length(remove)>0){
      theros = theros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  for (ii in 1:nrow(ichnotaxa)){
    remove = which(grepl(ichnotaxa[ii,],theros$gnl))
    if (length(remove)>0){
      theros = theros[-remove,]; # removing these occurrences
      removedoccs=removedoccs + length(remove)
    }
  }
  
  
  # There are some occurrences of Dinosauria not in any subclades, what are they?
  occids = setdiff(dinos$oid,c(sauros$oid,theros$oid,ornits$oid))
  dinoinx = NA;
  for (oo in 1:length(occids)){
    dinoinx[oo] = which(dinos$oid==occids[oo])
  }
  unique(dinos[dinoinx,]$mna)
  # A variety of species. The ones listed below are also removed from data.
  # Ichnotaxa not on Benson's list:
  # Harpedactylus (3 species,5 occs, gnn = 141796)
  # Yunnanpus     (1 species,1 occ , gnn = 91056)
  # Saurichnium   (4 species,3 occs, gnn = 141601)
  # Tetrapodium   (1 species,1 occ, gnn = 85805)
  #
  # Often referred to as nomen dubium, might be ichtyosaur.
  # Actiosaurus   (1 species,1 occ, gnn=268200)
  dinos = dinos[-which(dinos$gnn==141796),]
  dinos = dinos[-which(dinos$gnn==91056),]
  dinos = dinos[-which(dinos$gnn==141601),]
  dinos = dinos[-which(dinos$gnn==85805),]
  dinos = dinos[-which(dinos$gnn==268200),]
  
  # The data used for analysis is only from the Mesozoic, and these are
  # extracted inside the main analysis loop.

  ## == DONE CLEANING DATA ==
}


## -- Defining some variables and variable names --
Dinogroups <-c("Dinosauria","Ornithischia","Sauropodomorpha","Theropoda")
# Names for bins. PBDB intervals 112:138
interval.names =c("Maastrichtian","Campanian","Santonian","Coniacian","Turonian","Cenomanian","Albian","Aptian","Barremian","Hauterivian","Valanginian","Berriasian","Tithonian","Kimmeridgian","Oxfordian","Callovian","Bathonian","Bajacian","Aalenian","Toarcian","Pliensbachian","Sinemurian","Hettangian","Rhaetian","Norian","Carnian","Ladinian")
period.names = c("Triassic","Jurassic","Cretaceaous")
epoch.names  = c("Triassic","Early Jurassic","Mid Jurassic","Late Jurassic","Early Cretaceaous","Late Cretaceous")
# Defining stages.
# Interval years were collected from PBDB with manual scrutiny and correcting some obvious errors.
# Bins are [Start End Duration intervalindex (PBDB style), period index (1:3), epoch index (1:6)]
Bins = matrix(data=NA,nrow=27,ncol=6)
interval.names = array(data=NA,dim=c(27,1));
for (ss in 112:138){
  tmp=pbdb_interval(ss)
  Bins[ss-111,1]=tmp$eag
  Bins[ss-111,2]=tmp$lag
  interval.names[ss-111]=as.character(tmp$nam)
}
Bins[,3]=Bins[,1]-Bins[,2]# Durations
Bins[,4]=seq(112,138)# PBDB indices
Bins[,5] = c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1)# Periods: Triassic, Jurassic, Cretaceous
Bins[,6] = c(6,6,6,6,6,6,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,1,1,1,1)# Epoch: Triassic, Early Jurassic, Mid Jurassic, Late Jurassic, Early Cretaceaous, Late Cretaceous
# Note that Ladinian is here treated as part of Triassic as a whole.
rownames(Bins)<-interval.names
colnames(Bins)<-c("Start","End","Duration","PBDB_inx","Period_inx","Epoch_inx")
midpoints=(Bins[,1]+Bins[,2])/2

## ==== MAIN ANALYSIS LOOP ====
# norep is the number of replicated trials. Occurrences that span several intervals are
# stochastically placed in one of them with probabilities proportional to the length
# of the stages. All estimation are performed norep number of times and
# the median estimated sampling rate is used for further analysis.
if (doanalysis == TRUE){
  for (dd in 1:4){   # For all four datasets.
    # Switchin between datasets
    if (dd==1){
      usenow = dinos;
    } else if (dd==3){
      usenow = sauros;
    } else if (dd==2){
      usenow = ornits;
    } else if (dd==4){
      usenow = theros
    }
    
    # Occurrence counts. Temporary array to get number of unique species/genera
    if (dogenera==TRUE){
      tmp = createDataArrs_v4(usenow);
    } else {
      tmp = createDataArrs_v3(usenow);
    }
    
    # Generating the occurrence count arrays. The loop below generate norep number of matrices.
    tmpOccs = array(data=NA,dim=c(norep,nrow(tmp$Data),27)); # large array of all norep drawn occurrence matrices
    tmpSpec = matrix(data=NA,nrow=norep,ncol=27);  # large array of norep drawn occ.matrices' species richness counts.
    for (rr in 1:norep){
      # Replicating all collection of data.
      show(c(dd,rr))
      if (dogenera==TRUE){
        J = createDataArrs_v4(usenow); # v4 does the data by genera
      } else {
        J = createDataArrs_v3(usenow); # v3 does the data by species
      }
      tmpSpec[rr,]=colSums(J$Data>0); # Collecting species number
      tmpOccs[rr,,]=J$Data; # storing the occurrence count matrix
    }
    
    # Replicating sampling rate analysis for interval specific rates. This analysis
    # is performed with each replicated dataset.
    p_interval_now = array(data=NA,dim=c(norep,27,3));
    Spec_count_now = array(data=NA,dim=c(norep,27)); 
    for (rr in 1:norep){
      Occs = tmpOccs[rr,,]; # Now analysing this replicate [rr]
      Spec_count_now[rr,] = colSums(Occs>0); # Species count per interval
      for (ss in 1:27){
        if (sum((Occs[Occs[,ss]>0,ss])>1)){
          p_interval_now[rr,ss,]  <- estimatePoiss(array(Bins[ss,3],length(Occs[Occs[,ss]>0,ss])),Occs[Occs[,ss]>0,ss])
        }
      }
    }

    # Tallying range-through diversity. We assume that we have drawn enough replicates
    # to use the 'maximum' species number across replicates as range-through estimate. These will converge.
    tmp_RT = array(NA,dim=c(norep,27))
    for (rr in 1:norep){
      Occs = tmpOccs[rr,,];
      Range_through = array(NA,dim=dim(Occs))
      for (ss in 1:(dim(Occs)[1])){
        Range_through[ss,seq(min(which(Occs[ss,]>0)),max(which(Occs[ss,]>0)))]=1
      }
      tmp_RT[rr,] = colSums(Range_through,na.rm=TRUE); # Storing range-through richness count.
    }
    colnames(Occs)<-colnames(J$Data);
    rownames(Occs)<-rownames(J$Data);
    
    # Switchin between datasets, store to workspace
    if (dd==1){
      dinos_occs = Occs;
      dinos_occs_tmp = tmpOccs;
      dinos_times = J$Times;
      spec_count_dinos =  Spec_count_now;
      p_interval_dinos = p_interval_now;
      spec_count_rt_dinos = tmp_RT;
    } else if (dd==2){
      ornits_occs=Occs;
      ornits_occs_tmp = tmpOccs;
      ornits_times = J$Times;
      spec_count_ornits =  Spec_count_now;
      p_interval_ornits = p_interval_now;
      spec_count_rt_ornits = tmp_RT;
    } else if (dd==3){
      sauros_occs = Occs;
      sauros_occs_tmp = tmpOccs;
      sauros_times = J$Times;
      spec_count_sauros =  Spec_count_now;
      p_interval_sauros = p_interval_now;
      spec_count_rt_sauros = tmp_RT;
    } else if (dd==4){
      theros_occs=Occs;
      theros_occs_tmp = tmpOccs;
      theros_times = J$Times;
      spec_count_theros =  Spec_count_now;
      p_interval_theros = p_interval_now;
      spec_count_rt_theros = tmp_RT;
    }
    
  }
  
  
  ## Getting 'median' rates across replicates.
  p_int_median = array(NA,dim=c(27,4,3))
  for (ss in 1:27){
    # Dinosauria
    tix = which(p_interval_dinos[,ss,1]==sort(p_interval_dinos[,ss,1])[floor(length(sort(p_interval_dinos[,ss,1]))/2)])
    # if more than one hit, just use the first since they are identical
    p_int_median[ss,1,]=p_interval_dinos[tix[1],ss,]

    # Ornithischia
    tix = which(p_interval_ornits[,ss,1]==sort(p_interval_ornits[,ss,1])[floor(length(sort(p_interval_ornits[,ss,1]))/2)])
    # if more than one hit, just use the first since they are identical
    p_int_median[ss,2,]=p_interval_ornits[tix[1],ss,]
    
    # Sauropodomorpha
    tix = which(p_interval_sauros[,ss,1]==sort(p_interval_sauros[,ss,1])[floor(length(sort(p_interval_sauros[,ss,1]))/2)])
    # if more than one hit, just use the first since they are identical
    p_int_median[ss,3,]=p_interval_sauros[tix[1],ss,]
    
    # Theropoda
    tix = which(p_interval_theros[,ss,1]==sort(p_interval_theros[,ss,1])[floor(length(sort(p_interval_theros[,ss,1]))/2)])
    # if more than one hit, just use the first since they are identical
    p_int_median[ss,4,]=p_interval_theros[tix[1],ss,]
    
  }
  statnames = c('MLE','lower ci','upper ci')
  # Generating binomial sampling probability array
  dimnames(p_int_median)<-list(interval.names,Dinogroups,statnames)
  p_binos = array(NA,c(27,4,3));
  dimnames(p_binos)<-list(interval.names,Dinogroups,statnames)
  for (ss in 1:27){
    for (dd in 1:4){
      p_binos[ss,dd,]=1-exp(-p_int_median[ss,dd,]*Bins[ss,3])
    }
  }
  
  # Collecting observed species richness.
  SpecRich_test = array(NA,c(27,4));
  SpecRich_RT   = array(NA,c(27,4)); #Range through
  for (ss in 1:27){
    # Using the median species richness across replicate datasets.
    SpecRich_test[ss,1]=round(median(rowSums(dinos_occs_tmp[,,ss]>0)))
    SpecRich_test[ss,2]=round(median(rowSums(ornits_occs_tmp[,,ss]>0)))
    SpecRich_test[ss,3]=round(median(rowSums(sauros_occs_tmp[,,ss]>0)))
    SpecRich_test[ss,4]=round(median(rowSums(theros_occs_tmp[,,ss]>0)))
    SpecRich_RT[ss,1] = max(spec_count_rt_dinos[,ss])
    SpecRich_RT[ss,2] = max(spec_count_rt_ornits[,ss])
    SpecRich_RT[ss,3] = max(spec_count_rt_sauros[,ss])
    SpecRich_RT[ss,4] = max(spec_count_rt_theros[,ss])
  }
  dimnames(SpecRich_test)<-list(interval.names,Dinogroups)
  dimnames(SpecRich_RT)<-list(interval.names,Dinogroups)
  
  
  # Estimating true richnesses. 
  N_est_int = array(NA,dim=c(27,4,3))
  for (dd in 1:4){
    for (ss in 1:27){
      N_est_int[ss,dd,1]=floor(SpecRich_test[ss,dd]/p_binos[ss,dd,1]);
      N_est_int[ss,dd,2]=max(estimatetrue(SpecRich_test[ss,dd],p_binos[ss,dd,2]));
      N_est_int[ss,dd,3]=min(estimatetrue(SpecRich_test[ss,dd],p_binos[ss,dd,3]));
      
    }
  }
} else {
  if (dogenera==TRUE){
    load("RevisedAnalysis_Genera.RData")
  } else {
    load("RevisedAnalysis_Species.RData")
  }
}


## Below are other statistics in the main paper.

## -- Comparing binomial sampling probabilities with collection count ---
# Collection count is using the raw, unclean download, see before cleaning for calculation.
# Correlation between detrended log10 collections vs binomial probs. 
y = log10(nocol); # log10 number of dinosaur bearing collections across Mesozoic (log10(DBC))
residcol <- residuals(lm(y~midpoints)) # Residuals of a linear regression of (log10(DBC)) vs time
cor.test(p_binos[,1,1],residcol,use="pairwise.complete.obs")
cor.test(p_binos[,2,1],residcol,use="pairwise.complete.obs")
cor.test(p_binos[,3,1],residcol,use="pairwise.complete.obs")
cor.test(p_binos[,4,1],residcol,use="pairwise.complete.obs")

# Uncomment for figure plots
# par(mfrow=c(2,2))
# for (pp in 1:4){
#   plot(p_binos[,pp,1],residcol,type="p")
# }

# Richness weighted total number of species. We calculate a 'weighted' mean
# binomial probability (weighted by estimated richness) and used the no. observed
# total with this binomial prob. 
bino_tot = array(NA,c(4,3))
pois_tot = array(NA,c(4,3))
N_tot    = array(NA,c(4,3))
for (dd in 1:4){
  bino_tot[dd,]=c(sum(p_binos[,dd,1]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE),
                  sum(p_binos[,dd,2]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE),
                  sum(p_binos[,dd,3]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE))
  pois_tot[dd,]=c(sum(p_int_median[,dd,1]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE),
                  sum(p_int_median[,dd,2]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE),
                  sum(p_int_median[,dd,3]*N_est_int[,dd,1],na.rm=TRUE)/sum(N_est_int[,dd,1],na.rm=TRUE))
  
  
}
N_tot[1,] = c(estimatetrue(dim(dinos_occs_tmp)[2],bino_tot[1,1])[1], min(estimatetrue(dim(dinos_occs_tmp)[2],bino_tot[1,3])), max(estimatetrue(dim(dinos_occs_tmp)[2],bino_tot[1,2])))
N_tot[2,] = c(estimatetrue(dim(ornits_occs_tmp)[2],bino_tot[2,1])[1], min(estimatetrue(dim(ornits_occs_tmp)[2],bino_tot[2,3])), max(estimatetrue(dim(ornits_occs_tmp)[2],bino_tot[2,2])))
N_tot[3,] = c(estimatetrue(dim(sauros_occs_tmp)[2],bino_tot[3,1])[1], min(estimatetrue(dim(sauros_occs_tmp)[2],bino_tot[3,3])), max(estimatetrue(dim(sauros_occs_tmp)[2],bino_tot[3,2])))
N_tot[4,] = c(estimatetrue(dim(theros_occs_tmp)[2],bino_tot[4,1])[1], min(estimatetrue(dim(theros_occs_tmp)[2],bino_tot[4,3])), max(estimatetrue(dim(theros_occs_tmp)[2],bino_tot[4,2])))

## ===  Printing figures === ##
# Uncomment three lines above and one below the "source('xxx.R')" to have them print to pdf.
if (dofigures==TRUE){
  dospec = sum(dogenera==FALSE); #yaxis label 1- species richness 2 - genus richness

  #  -- Sampling rates (figure 1 in paper) -- 
  docis=FALSE # if plotting the confidence interval of the sampling data.
  # filname = paste('Fig_1_Sampling_Species_',as.character(Sys.Date()),'.pdf',sep="")
  # if (dogenera==TRUE){filname = gsub(pattern='Species',replacement='Genera',filname)}
  # pdf(file=filname,pointsize=42,height=12,width=24)
  source('Figure_1_Sampling_191115.R')
  # dev.off()
  
  
  # -- Sampling rates with confidence intervale (in ESM) --
  docis=TRUE # Switch for printing with confidence intervals
  # filname = paste('Fig_SI4_Sampling_Species_CI_',as.character(Sys.Date()),'.pdf',sep="")
  # if (dogenera==TRUE){filname = gsub(pattern='Species',replacement='Genera',filname)}
  # png(file="Fig_ESM_samplingSpec.png",pointsize=42,height=12,width=24,units="in",res=300)
  source('Figure_1_Sampling_191115.R')
  # dev.off()

  # -- Estimated richness (figure 2 in main paper) --
  # filname = paste('Fig_3_SpeciesRichness_',as.character(Sys.Date()),'.pdf',sep="")
  # if (dogenera==TRUE){filname = gsub(pattern='Species',replacement='Genera',filname)}
  # png(file="Genusrichness.png",pointsize=42,height=24,width=32,units="in",res=600)
  source('RichnessFig_rev.R')
  # dev.off()
  
  
  # -- Sampling rates over the replicated 100 datasets --
  # filname = paste('Fig_SI5_ReplicateSampling_',as.character(Sys.Date()),'.pdf',sep="")
  # if (dogenera==TRUE){filname = gsub(pattern='Species',replacement='Genera',filname)}
  # png(file='Fig_ESM_Replicatedsampling.png',pointsize=32,width=16,height=12,units="in",res=600)
  # pdf(file=filname,pointsize=32,width=24,height=24)
  source('plotallsamplingrates_SI.R')
  # dev.off()
  
}