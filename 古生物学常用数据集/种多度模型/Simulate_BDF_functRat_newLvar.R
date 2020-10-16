SimulateBDF <- function(spec,mu,lambda,lambdavar,tmax,n_init,dt){
  ## Simulating a birth-death-fossilize model.
  # Written by Jostein Starrfelt (jostein.starrfelt[AT]ibv.uio.no)
  # Last checked 09.09.15
  # Should be updated to only use 'functional' approach to rates actually defined as function. Now
  # if any are defined as functions, all are treated as functions. This is inefficient. [JOS 051015]
  # Also there must be an error here;
  
  #   
  # Updated to have potentially variable rates over time. These rates are inputted as functions taking the argument t-time.
  # This now works but then without lambdavar. Implementing lambdavar as a function that outputs a scaling factor for of sampling. 
  # Lambdavar is then essentially a scaling distribution, which should have mean 1. The sampling rate for the whole clade is
  # lambda*lambdavar.
  # v3:
  # In this way one can implement whatever kind of distribution of differential sampling rates across lineages. The sampling rate for an
  # individual lineage at time t will then be lambda(t) * lambdavar(1). Future implementations might include also this as a function of time.
  # If not species do not vary in sampling rate, set lambdavar to 1. This is different from previous applications. Both the 'fixed' rate and
  # the functional code below needs to be augmented.
  
  # Implementing in fixed rate code; DONE.
  
  # Implementing in variable rate code; DONE, not tested.
  
  # Calling functions
  # When calling a function you can specify arguments by position, by complete name, or by partial name. Arguments are matched first by exact name (perfect matching), then by prefix matching, and finally by position. 
  
  # How to implement
  #   spec = 0.01;
  #   mu = 0;
  #   lambda = 0.4;
  #   lambdavar = 0;
  #   tmax = 10;
  #   n_init = 100;
  # Splitting code into fixed and variable rates
  # Assuming fixed speciation and extinction rates for whole interval.
  if (any(list(typeof(spec),typeof(mu),typeof(lambda)) == "closure")){
    # Then one or more of these arguments are functions and needs to be treated accordingly.
    # If one or more are given as constants; make into functions;
    if (typeof(spec)=="double"){spold = spec; spec <- function(t){spold}; spec <- Vectorize(spec)}
    if (typeof(mu)=="double"){muold = mu; mu <- function(t){muold}; mu <- Vectorize(mu)}
    if (typeof(lambda)=="double"){lold = lambda; lambda <- function(t){lold}; lambda <- Vectorize(lambda)}
    # The following lines calculates the estimated mean number of cumulative species richness according to e BD process;
    # cf eq 41 in Kendall 1948
    myrhotmp <-  function(t){(mu(t)-spec(t))};
    myrho    <- function(t) {integrate(myrhotmp,0,t)}
    tmpfun <- function(t){exp(-myrho(t)$value)}
    mytmp  <- function(t){exp(-myrho(t)$value*spec(t))}
    nosp <- round(n_init*(1 + integrate(mytmp,0,tmax)$value))
    Out = array(NA,c(nosp*100,3));
    # Should put in dynamic checker for this array too, since it might be too small sometimes.
    Lamds = array(NA,c(nosp*100,1)); # for testing the lambdavalues
    done = 0; #checker for done calculating or not.
    Out[1:n_init,1] = 0; # initial lineages have start at t=0;
    tix = 1; # current lineage being simulated
    ntix = n_init+1; # current index into first entry in Out not populated with species yet.
    while (done==0){
      # while not done.
      # Here I will insert a switch to redefine the functions spec, mu and lambda, in cases where they will vary between individual lineages.
      # Exactly how this will be simulated remains to be seen.
      # if (doindfun==1){
      #redefine functions so they are lineage specific here. Most relevant for sampling function
      # }
      # Now simulating lineage tix.
      dead = 0; # switches to 1 when this lineage goes extinct;
      t_x = Out[tix+1-1,1]; # simulated time, this starts at time of speciation for lineage tix
      # dt = 1e-2; # Traveling along time with this resolution. The code below will integrate the relevant rates over many short intervals
      # of this duration.
      while (dead==0){
        if (runif(1)<integrate(mu,t_x,t_x+dt)$value){
          # if this is true it has died
          dead = t_x;
          Out[tix,2] = t_x;
          
        } else{
          t_x = t_x+dt; # just update time
        }
        if (t_x>tmax){
          # If not dead by end of simulation it has not died, but has duration until time tmax
          dead=-1; # To stop loop
          Out[tix,2] = tmax;  # This lineage has not gone extinct before tmax.
        }
        
      }
      # Now we have the full duration for this lineage.
      # Getting the speciation event for this lineage;
      # List of speciation times;
      S = array(NA,10*round(integrate(spec,0,tmax)$value));
      stix = 1; # Entry into the S array
      # dt = 1e-3;
      t_x = Out[tix+2-2,1]; # WHen this lineage first appears
      while (t_x<=Out[tix,2]){
        if (runif(1)<integrate(spec,t_x,t_x+dt)$value){
          # if this is true a speciation events happens now.
          # Need to check if S is big enough;
          if (stix>length(S)){
            S_tmp = S;
            S = array(NA,dim=c(round(length(S_tmp)*1.2),1)); # increasing size by 20 %
            S[1:stix-1] = S_tmp[1:stix-1];
          }
          S[stix] = t_x; # speciation time
          stix = stix+1; # increase counter
          t_x = t_x+dt;  # increase time
        } else {
          # nothing happens, time is updated
          t_x = t_x + dt;
        }
      }
      S = S[!is.na(S)]; # keeping only true speciation events from preallocated array.
      # show(S)
      # show(tix)
      
      if ((ntix+length(S)-1)>nrow(Out)){
        ## need to enlarge OUT
        oldout= Out;
        # Generating new array ~50 % bigger.
        Out = array(NA,c(round(nrow(Out)*1.5),3));
        Out[1:(ntix-1),]=oldout[1:(ntix-1),];
      }
      
      if (length(S)>0){
        Out[ntix:(ntix+length(S)-1)]=S;
        ntix =ntix +length(S); # updating index into first non-populated part of Out.
      }
      
      
      
      # Fossilization/Sampling; 
      nofos = 0; # counting number of fossils for this lineage
      t_x = Out[tix+3-3,1]; #starting at first appearance of this lineage
      # Lambda scaling factor for this lineage
      lambdatmp = lambdavar(1);
      Lamds[tix] = lambdatmp;
      while (t_x<=Out[tix,2]){
        # show(t_x)
        # For the duration of this lineage
        if (runif(1)<lambdatmp*integrate(lambda,t_x,t_x+dt)$value){
          # If fossilization event inside the time t_x ... t_x + dt, adda fossil to count
          nofos = nofos+1;
          t_x = t_x + dt; # update time
        } else {
          # No fossilization event
          t_x = t_x + dt; # update time.
        }
      }
      Out[tix,3] = nofos;
      
      tix = tix+1; #now done simulating lineage tix, moving to tix+1
      # show(tix/ntix)
      if (tix==(ntix)){
        # if this was the last lineage in the array, we are done.
        done=1;
      }
    }
    
    
  } else {
    # If rates are fixed over time.
    rho <- function(t,spec,mu){mu*t - spec*t}
    rho2 <- function(t){mu*t - spec*t}
    tmpfun <- function(t){exp(-rho2(t))*spec}
    nosp <- function(t,spec,mu,n_init){ n_init*(1 + integrate(tmpfun,0,t)$value)}
    # Rough estimate of total species richness expected by Kendall
    minspec <- round(nosp(tmax,spec,mu,n_init))
    Out = array(NA,c(minspec*10,3)); # preallocating the array for species
    Lamds = array(NA,c(minspec*10,1)); # for testing the lambdavalues
    done = 0;
    Out[1:n_init,1]=0;
    tix = 1;
    ntix = n_init+1; #index into first row in Out with no entries yet.
    while (done==0) {
      # When does this lineage go extinct
      if (mu>0){ # if extinction rate is nonzero
        Out[tix,2] = min(Out[tix,1]+rexp(1,mu),tmax)
      } else {
        Out[tix,2] = tmax;
      }
      # Drawing number of fossils for this lineage
      # draw one number from Poisson distribution with mean drawn from a normal distribution with mean lambda
      # and st.dev lambdavar times the duration of the taxon.
      # Out[tix,3] = rpois(1,max(0,rnorm(1,lambda,lambdavar))*(Out[tix,2]-Out[tix,1]))
      # This was before lambdavar = 1 if no variability.
      indsamp = lambda*lambdavar(1);
      Out[tix,3] = rpois(1,max(0,indsamp*(Out[tix,2]-Out[tix,1]))); #
      Lamds[tix] = indsamp/lambda; # Now this is also a scaling factor
      # Drawing waiting times to possible speciation events.
      if (spec>0){
        tmp <- rexp(1e1,spec)
        tmptix = 1; 
        # The above draw might be too small;
        while (sum(tmp)<(Out[tix,2]-Out[tix,1])){
          tmp = rexp(10^(tmptix),spec);
          tmptix = tmptix+1; # increase number of draws until end of
          # lineage is traversed
        }
        spectimes = cumsum(tmp)<(Out[tix,2]-Out[tix,1]);
        if (sum(cumsum(tmp)<(Out[tix,2]-Out[tix,1]))>0){
          
          if ((ntix+sum(spectimes)-1)>nrow(Out)){
            ## need to enlarge OUT
            oldout= Out;
            # Generating new array ~50 % bigger.
            Out = array(NA,c(round(nrow(Out)*1.5),3));
            Out[1:(ntix-1),]=oldout[1:(ntix-1),];
          }
          # if some of these 'speciation times' are inside the actual duration of taxon tix
          Out[ntix:(ntix+sum(spectimes)-1),1]=Out[tix,1]+cumsum(tmp)[spectimes]
          ntix = ntix+sum(spectimes)
          
        }
      } else {
        # tmp = 0;
        # Do nothing if speciation rate is 0
      }
      
      tix = tix+1;
      if (tix==(ntix)){
        done=1;
      }
      
    }
  }
  Out = Out[!is.na(Out[,1]),]; # removin excess NA's
  Lamds = Lamds[!is.na(Lamds)];
  return(list(Out,Lamds))
}
