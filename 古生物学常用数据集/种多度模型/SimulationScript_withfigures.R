# Script to perform simulations used in the TRiPS main paper.

rm(list=ls()) # emptying workspace
dosimulation = TRUE # If set to true new simulations are run. If false, figures for paper are plotted.
# If rerunning simulations, set limits to parameters below.


load('RevisedAnalysis_Species.Rdata') # Loading in estimation results for figure 3 (main paper).
if (dosimulation==FALSE){
  load('Simus_run_220915.RData')
  
  
} else {
  source('Simulate_BDF_functRat_newLvar.R')
  library(stats4)
  source('functions_DinoBino_toweb.R')
  library(lhs)
  norep = 1e3; #Number of simulations.
  set.seed(10) #Setting seed for reproducability
  Parset = array(NA,c(6,3))
  # Minimum, maximum and type of distribution for the 5 parameters;
  # Speciation rate :  spec
  Parset[1,1] = 1e-3; # minimum
  Parset[1,2] = 2e-3; # maximum
  Parset[1,3] = 2;    # 1 for uniform 2 for log-uniform
  # Exctinction rate :  mu
  Parset[2,1] = 1e-3; # minimum
  Parset[2,2] = 2e-3; # maximum
  Parset[2,3] = 2;    # 1 for uniform 2 for log-uniform
  # Sampling rate : lambda
  Parset[3,1] = 1e-3;# minimum
  Parset[3,2] = 0.5; # maximum
  Parset[3,3] = 1;   # 1 for uniform 2 for log-uniform
  # Individual variability in sampling rate :  varlam
  Parset[4,1] = 0;  # minimum
  Parset[4,2] = .3; # maximum
  Parset[4,3] = 1;  # 1 for uniform 2 for log-uniform
  # Number of initial lineages
  Parset[6,1] = 10;  # minimum
  Parset[6,2] = 250; # maximum
  Parset[6,3] = 1;   # 1 for uniform 2 for log-uniform
  # Duration of simulation
  Parset[5,1] = 2;  # minimum
  Parset[5,2] = 20; # maximum
  Parset[5,3] = 1;  # 1 for uniform 2 for log-uniform
  
  # Drawing parameters latinsquare hypercube style.
  tmp = randomLHS(norep,6)
  Pars = array(NA,c(norep,6));
  for (ii in 1:6){
    if (Parset[ii,3]==1){ #uniform
      Pars[,ii] = Parset[ii,1] + (Parset[ii,2]-Parset[ii,1])*tmp[,ii]
    } else { # Assume log uniform
      Pars[,ii] = 10^(log10(Parset[ii,1]) + (log10(Parset[ii,2])-log10(Parset[ii,1]))*tmp[,ii])
    }
  }
  Pars[,6]=round(Pars[,6]); #rounding number of initial lineages
  
  par(mfrow=c(3,2))
  for (ii in 1:6){
    hist(Pars[,ii])
  }
  Occs_Big = array(0,c(norep,1e3)); # Assume there are no species with > 1000 occurrences. All above 1000 are put in last bin.
  N_big    = array(NA,c(norep,6));
  L_big    = array(NA,c(norep,4));
  for (rr in 1:norep){
    show(rr/norep)
    lambdavar = function(n) runif(1,1-Pars[rr,4],1+Pars[rr,4])
    tmax = Pars[rr,5];
    # pre 210215
    J = SimulateBDF(Pars[rr,1],Pars[rr,2],Pars[rr,3],lambdavar,Pars[rr,5],Pars[rr,6])
    # J2 = SimulateBDF(Pars[rr,1],Pars[rr,2],Pars[rr,3],lambdavar,Pars[rr,5],Pars[rr,6])
    # <- function(spec,mu,lambda,lambdavar,tmax,n_init){
      
    
    Out = J[[1]]
    N_big[rr,1] = length(Out[,1]); # true richness
    N_big[rr,2]  = sum(Out[,3]>0); # number of observed species
    tmp <- get_ltt(Out);
    N_big[rr,3]= max(tmp[[2]]); # maximum number of lineages at the same time.
    if (doanalysis == TRUE){
      use = Out[,3]>0;
      TRiPS = array(NA,c(1,3))
      if (any(Out[use,3]>1)){
        p_pois = NA;
        p_pois = tryCatch(estimatePoiss(rep(tmax,sum(use)),Out[use,3])) #Saving estimated rates
        if (!is.na(p_pois)){ #then it worked
          p_bino = 1-exp(-p_pois*tmax) #making them into binomial probs
          TRiPS[1] =  floor(sum(Out[,3]>0)/p_bino[1]);
          TRiPS[2] = min(estimatetrue(sum(Out[,3]>0),p_bino[3]))
          TRiPS[3] = max(estimatetrue(sum(Out[,3]>0),p_bino[2]))
        }
      }
      N_big[rr,4:6] = TRiPS;
      L_big[rr,1:3] = p_pois;
      L_big[rr,4]  = mean(J[[2]])*Pars[rr,3]; # mean sampling rate among lineages. NOTE THIS IS NOT USEFUL FOR FUNCTIONAL SAMPLING.
      tmp = rle(sort(Out[Out[,3]>0,3]))
      data = array(NA,c(length(tmp$lengths),2))
      data[,1] = tmp$values
      data[,2] = tmp$lengths
      Occs_Big[rr,data[,1]] = data[,2]; # storing the occurrence counts/freq of freqs
    }
  }
  
  use = !is.na(N_big[,4])
  # Quick calc:
  # Success rate
  # Correlation
  # Mean scaled error
  sum(N_big[use,1]<=N_big[use,6] & N_big[use,1]>=N_big[use,5])/sum(use)
  cor(N_big[use,1],N_big[use,4],use="pairwise.complete.obs")
  sum((N_big[use,1]-N_big[use,4])/(N_big[use,1]))/length(use)
  
  # Doing calculations for making figure 3 in main paper.
  xs = seq(0,Parset[3,2],by=0.05)
  ys = seq(2,Parset[5,2],by=1)
  hits = array(NA,c(length(xs)-1,length(ys)-1))
  hitslam = array(NA,c(length(xs)-1,length(ys)-1))
  cors = array(NA,c(length(xs)-1,length(ys)-1))
  corslog = array(NA,c(length(xs)-1,length(ys)-1))
  cors_obs = array(NA,c(length(xs)-1,length(ys)-1))
  SME  = array(NA,c(length(xs)-1,length(ys)-1)) # bias
  SRMSE  = array(NA,c(length(xs)-1,length(ys)-1))  # scaled root mean square error
  nos = array(NA,c(length(xs)-1,length(ys)-1))  # scaled root mean square error
  
  for (ii in 1:(length(xs)-1)){
    for (jj in 1:(length(ys)-1)){
      
      use = Pars[,3]<xs[ii+1] & Pars[,3]>xs[ii] & Pars[,5]<ys[jj+1] & Pars[,5]>=ys[jj] & !is.na(N_big[,4])
      if (sum(use)>1){
      nos[ii,jj] = sum(use)
      hits[ii,jj] = sum(N_big[use,1]<=N_big[use,6] & N_big[use,1]>=N_big[use,5],na.rm=TRUE)/sum(use)
      hitslam[ii,jj] = sum(L_big[use,4]<=L_big[use,3] & L_big[use,4]>=L_big[use,2])/sum(use)
      cors[ii,jj] = cor(N_big[use,1],N_big[use,4],use="pairwise.complete.obs")
      corslog[ii,jj] = cor(log10(N_big[use,1]),log10(N_big[use,4]),use="pairwise.complete.obs")
      cors_obs[ii,jj] = cor(N_big[use,1],N_big[use,2]);
      SME[ii,jj]  = sum((N_big[use,1]-N_big[use,4])/(N_big[use,1]))/sum(use)
      SRMSE[ii,jj] = sum(((N_big[use,1]-N_big[use,4])^2)/(N_big[use,1]^2))/sum(use)^2}
    }
  }
  
  xpl = (xs[1:(length(xs)-1)] + xs[2:length(xs)])/2
  ypl = (ys[1:(length(ys)-1)] + ys[2:length(ys)])/2
}

# Printing figure 3 in main paper.
## 

tcex = 0.8
# printfilename = 'Figure3_rev.pdf'
# pdf(printfilename,width=16,height=16,pointsize=24)
tmpcols = grey.colors(11)
layout(matrix(c(0,0,1,2),nrow=2,ncol=2,byrow=TRUE),c(5,1),c(1,5))
plot(p_int_median[,1,1],Bins[,3],pch=15,col=color.list[ii],xlim = c(Parset[3,1]-2e-3,Parset[3,2]),ylim = c(Parset[5,1]*0.98,Parset[5,2]),xaxs="i",yaxs="i",xlab="",ylab="")
# xaxs = "i" forces the +/- 6% addition to the limits to be ignored, i.e. limits are ACTUAL limits in the graph
for (ii in 1:(length(xs)-1)){
  for (jj in 1:(length(ys)-1)){
    rect(xs[ii],ys[jj],xs[ii+1],ys[jj+1],col=tmpcols[round(10*hits[ii,jj]+1)],border=NA)
  }
}
title( xlab=expression(paste(lambda," - sampling rate ")), ylab='Interval duration (Ma)',cex.lab=tcex)
for (ii in 1:4){
  points(p_int_median[,ii,1],Bins[,3],pch=(14+ii),col=color.list[ii])
}
xbar = seq(0,1,by=0.1)
plot(0,0,xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
for (ii in 1:(length(xbar)-1)){
  rect(0,xbar[ii],1,xbar[ii+1],col=tmpcols[ceiling(10*xbar[ii])+1],xaxt="n")
}
axis(2,las=2)
# dev.off()


## NOTE! THe L_big has entries where it should not have;
# For cases where TRiPS fail, the code has for some reason stored last rounds lambda estimate!
tmp=which(is.na(N_big[,4]))
L_big[tmp,1]-L_big[tmp-1,1]

L_big[tmp,]=NA
nix = 21;
bigX = array(NA,c(nix,6)); # x values for plots for each parameter
for (pp in 1:6){
  if (Parset[pp,3]==1){
    bigX[,pp]=seq(Parset[pp,1],Parset[pp,2],length=nix)
  } else {
    bigX[,pp] = 10^(seq(log10(Parset[pp,1]),log10(Parset[pp,2]),length=nix))
  }
}

# And we want to calculate the fraction of 'hits', the correlation and the mean scaled error
# across these 6 parameters.
Output = array(NA,c(6,nix,4)); # Parameter by x by [hits,correlation,mse]
Nouse  = array(NA,c(6,nix));   # number of simulations in this bin
Output_L = array(NA,c(6,nix,4)); # Parameter by x by [hits,correlation,mse]
Nouse_L  = array(NA,c(6,nix));   # number of simulations in this bin
Output_B = array(NA,c(6,nix,4)); # Parameter by x by [hits,correlation,mse]
Nouse_B  = array(NA,c(6,nix));   # number of simulations in this bin
rownames(Parset)<-c("Speciation rate","Extinction rate","Sampling rate","Sampling variability","Duration","No initial lineages")
# Prover mean (unscaled) error for binomials, since they are bounded
Groups = array(NA,c(100000,6));
for (pp in 1:6){
  for (ix in 1:(nix-1)){
    
    use = which(Pars[,pp]>=bigX[ix,pp] & Pars[,pp]<bigX[ix+1,pp] & !is.na(N_big[,4]))
    Groups[use,pp]=ix;
    Nouse[pp,ix]=length(use);
    
    hit = N_big[use,5]<=N_big[use,1] & N_big[use,6]>=N_big[use,1]
    Output[pp,ix,1] = sum(hit)/length(use);
    Output[pp,ix,2] = cor(log10(N_big[use,1]),log10(N_big[use,4]),use="pairwise.complete.obs")
    Output[pp,ix,4] = cor(log10(N_big[use,1]),log10(N_big[use,2]),use="pairwise.complete.obs")
    Output[pp,ix,3] = sum((N_big[use,1]-N_big[use,4])/N_big[use,1])/length(use)
    
    hitL = L_big[use,2]<=Pars[use,3] & L_big[use,3]>=Pars[use,3]
    Output_L[pp,ix,1] = sum(hitL)/length(use);
    Output_L[pp,ix,2] = cor(Pars[use,3],L_big[use,1],use="pairwise.complete.obs")
    Output_L[pp,ix,3] = sum((Pars[use,3]-L_big[use,1])/Pars[use,3])/length(use)    
    tmpestbinl = 1-exp(-L_big[use,2]*Pars[use,5]);
    tmpestbinh = 1-exp(-L_big[use,3]*Pars[use,5]);
    tmpestbin = 1-exp(-L_big[use,1]*Pars[use,5]);
    truebin = N_big[use,2]/N_big[use,1];
    hitB = (tmpestbinl<truebin & tmpestbinh >truebin);
    Output_B[pp,ix,1] = sum(hitB)/length(use)
    Output_B[pp,ix,2] = cor(truebin,tmpestbin,use="pairwise.complete.obs")
    Output_B[pp,ix,3] = sum((truebin-tmpestbin))/length(use)    
  }
}

# Plotting success rate, correlation and mean scaled error of RICHNESS estimates vs parameter ranges.
# png('SI_fig_1.png',width=16,height=16,units="cm",pointsize=12,res=600)
par(mfrow=c(3,2))
par(mar = c(5,5,2,5))
for (pp in 1:6){
  if (Parset[pp,3]==2){
    plot(bigX[,pp],Output[pp,,1],type="p",pch=15,log="x",ylim=c(0,1),ylab="Correlation, success rate",xlab=rownames(Parset)[pp])
    lines(bigX[,pp],Output[pp,,2],type="p",pch=16,col="green")
    title(main="Richness")
    # lines(bigX[,pp],Output[pp,,4],type="l",col="blue")
    
    par(new=TRUE)
    plot(bigX[,pp],Output[pp,,3],log="x",type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.05,0.2))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  } else {
    plot(bigX[,pp],Output[pp,,1],type="p",pch=15,ylim=c(0,1),xlab=rownames(Parset)[pp],ylab="Correlation, success rate")
    lines(bigX[,pp],Output[pp,,2],type="p",pch=16,col="green")
    title(main="Richness")
    # lines(bigX[,pp],Output[pp,,4],type="l",col="blue")
    
    par(new=TRUE)
    plot(bigX[,pp],Output[pp,,3],type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.05,0.2))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  }
  
  
}
# dev.off()
# 
# Plotting success rate, correlation and mean scaled error of SAMPLING RATE estimates vs parameter ranges.
# png('SI_fig_2.png',width=16,height=16,units="cm",pointsize=12,res=600)
par(mfrow=c(3,2))
par(mar = c(5,5,2,5))
for (pp in 1:6){
  if (Parset[pp,3]==2){
    plot(bigX[,pp],Output_L[pp,,1],type="p",pch=15,log="x",ylim=c(0,1),ylab="Correlation, success rate",xlab=rownames(Parset)[pp])
    lines(bigX[,pp],Output_L[pp,,2],type="p",pch=16,col="green")
    title(main="Sampling rate")
    # lines(bigX[,pp],Output[pp,,4],type="l",col="blue")
    
    par(new=TRUE)
    plot(bigX[,pp],Output_L[pp,,3],log="x",type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.5,0.5))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  } else {
    plot(bigX[,pp],Output_L[pp,,1],type="p",pch=15,ylim=c(0,1),xlab=rownames(Parset)[pp],ylab="Correlation, success rate")
    lines(bigX[,pp],Output_L[pp,,2],type="p",pch=16,col="green")
    title(main="Sampling rate")
    # lines(bigX[,pp],Output[pp,,4],type="l",col="blue")
    
    par(new=TRUE)
    plot(bigX[,pp],Output_L[pp,,3],type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.5,0.5))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  }
  
  
}
# dev.off()

# Plotting success rate, correlation and mean scaled error of SAMPLING PROBABILITY estimates vs parameter ranges.
par(mfrow=c(3,2))
par(mar = c(5,5,2,5))
for (pp in 1:6){
  if (Parset[pp,3]==2){
    plot(bigX[,pp],Output_B[pp,,1],type="p",pch=15,log="x",ylim=c(0,1),ylab="Correlation, success rate",xlab=rownames(Parset)[pp])
    lines(bigX[,pp],Output_B[pp,,2],type="p",pch=16,col="green")
    title(main="Sampling prob")
    par(new=TRUE)
    plot(bigX[,pp],Output_B[pp,,3],log="x",type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.5,0.5))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  } else {
    plot(bigX[,pp],Output_B[pp,,1],type="p",pch=15,ylim=c(0,1),xlab=rownames(Parset)[pp],ylab="Correlation, success rate")
    lines(bigX[,pp],Output_B[pp,,2],type="p",pch=16,col="green")
    title(main="Sampling prob")
    par(new=TRUE)
    plot(bigX[,pp],Output_B[pp,,3],type="p",pch=17,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-.5,0.5))
    axis(4,col="red")
    mtext("Mean Scaled Error",side=4,line=3,col="red",cex=0.8)
  }
  
  
}
# dev.off()
# 
# 
