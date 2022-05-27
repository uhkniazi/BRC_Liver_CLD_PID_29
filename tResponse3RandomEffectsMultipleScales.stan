data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1> Nclusters3; 
  int<lower=1> Nclusters4;
  int<lower=1> NScaleBatches1; // number of batches of scale terms for group 1 (subgroups in group 1) 
  int<lower=1> NScaleBatches2; // number of batches of scale terms for group 2 (subgroups in group 2) 
  int<lower=1> NScaleBatches3;
  //int<lower=1> NScaleBatches4;
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each group 1 parameter to a data point 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each group 2 parameter to a data point 
  int<lower=1, upper=Nclusters3> NgroupMap3[Ntotal]; 
  int<lower=1, upper=Nclusters4> NgroupMap4[Ntotal]; 
  int<lower=1, upper=NScaleBatches1> NBatchMap1[Nclusters1]; // expanding vector to map each scale term to the relevant 
                                                    // set of coefficients from group 1
  int<lower=1, upper=NScaleBatches2> NBatchMap2[Nclusters2]; // expanding vector to map each scale term to the relevant 
                                                    // set of coefficients from group 2                                                    
  int<lower=1, upper=NScaleBatches3> NBatchMap3[Nclusters3];
  //int<lower=1, upper=NScaleBatches4> NBatchMap4[Nclusters4];
  real y[Ntotal]; // response variable
  vector[Ntotal] X; // slope variable
  // additional parameters
  real intercept;
  real intercept_sd;
  int<lower=1> Nnu; // number of nu terms for each subset of observations 
  int<lower=1> NsigmaPop;
  int<lower=1, upper=Nnu> NnuMap[Ntotal]; // mapping variable 
  int<lower=1, upper=NsigmaPop> NsigmaPopMap[Ntotal]; // mapping variable 
}

// transformed data {
  // }

parameters {
  // parameters to estimate in the model
  real betas; // constant intercept term
  real slope; // population slope term
  vector<lower=0.01>[NScaleBatches1] sigmaRan1; // random effect standard deviations for sub-batches in group 1
  vector<lower=0.01>[NScaleBatches2] sigmaRan2; // random effect standard deviations for sub-batches in group 2
  vector<lower=0.01>[NScaleBatches3] sigmaRan3; 
  //vector<lower=0.01>[NScaleBatches4] sigmaRanSlope;  
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
  vector[Nclusters3] rGroupsJitter3;
  vector[Nclusters4] rGroupsSlopes;
  real<lower=1> nu[Nnu]; // normality parameter for t distribution or degree of freedom 
  real<lower=0.01> sigmaPop[NsigmaPop]; // population standard deviation
}

transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] newSlope;
  newSlope = slope + rGroupsSlopes[NgroupMap4];
  // fitted value
  for (i in 1:Ntotal){
    mu[i] = X[i] * newSlope[i];
  }
  // fitted values
  mu = mu + betas + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2] + rGroupsJitter3[NgroupMap3];
}

model {
  vector[Nclusters1] sigmaRan1_expanded;
  vector[Nclusters2] sigmaRan2_expanded;
  vector[Nclusters3] sigmaRan3_expanded;
  // hyperparameters
  sigmaRan1 ~ exponential(1);
  sigmaRan2 ~ exponential(1);
  sigmaRan3 ~ exponential(1);
  sigmaPop ~ exponential(1);
  nu ~ exponential(1/29.0);
  
  // regression coefficients
  betas ~ normal(intercept, intercept_sd);
  slope ~ normal(0, 2);
  // parameter mapping to longer array, as each batch of coefficients map separately to a unique scale term
  sigmaRan1_expanded = sigmaRan1[NBatchMap1];
  sigmaRan2_expanded = sigmaRan2[NBatchMap2];
  sigmaRan3_expanded = sigmaRan3[NBatchMap3];
  // varying coefficients centered at 0
  rGroupsJitter1 ~ normal(0, sigmaRan1_expanded);
  rGroupsJitter2 ~ normal(0, sigmaRan2_expanded);
  rGroupsJitter3 ~ normal(0, sigmaRan3_expanded);
  rGroupsSlopes ~ normal(0, 2);
  // likelihood function
  y ~ student_t(nu[NnuMap], mu, sigmaPop[NsigmaPopMap]);  
}
