# Write time series matrix of cmip6 annual averages.
# Data is initially processed by Andy Wiltshire.
# dougmcneall@gmail.com

# In this version of the code, the weighting scheme doesn't affect the model projections
# a great deal. It's not clear if that's because i've messed up coding the weighting scheme,
# if there isn't much information in the temperature scheme, or if the models are all equally
# dissimilar from the observations.
#
# Some tactics for testing this might be:
# Test the weighting scheme using each model as an observation.
# Try using different free parameters in the weighting scheme.


library(ncdf4)

# The tricky part is creating a list of all the relevant files
# in the directory.

# r = realization index
# i = initialization index
# p = physics index
#f = forcing index

extractModelnames = function(filenames, startstring, stopstring){
  # Function to extract model names from a vector of cmip filenames
  # that have similar startstring and stopstring

  fs = lapply(filenames,regexpr, pattern = startstring)
  fe = lapply(filenames,regexpr, pattern = stopstring)
  out = rep(NA, length(filenames))

  for(i in 1:length(filenames)){
    out[i] = substr(filenames[i], attr(fs[i][[1]], 'match.length')+1, fe[[i]][1]-1)
  }
  out
}

# -----------------------------------------------------------------------------------------------
# Use historical tas compared with hardcrut observations to weight the models.
# Projections are using SSP585 (for maximum signal).
# -----------------------------------------------------------------------------------------------

# data directories
fn.ssp585 = dir('/data/users/hadaw/cmip6/areaavg/tas/', pattern = 'ssp585_r1i1p1')
fn.historical = dir('/data/users/hadaw/cmip6/areaavg/tas/', pattern = 'historical_r1i1p1')

# model name lists
mods.ssp585 = extractModelnames(fn.ssp585, startstring = 'tas_Amon_', stopstring = "_ssp585")
mods.historical = extractModelnames(fn.historical, startstring = 'tas_Amon_',
  stopstring = "_historical")

# Find models common to historical and future projection
common.mods = intersect(mods.ssp585, mods.historical)

fn.ssp585.kept = fn.ssp585[match(common.mods, mods.ssp585)]
fn.historical.kept = fn.historical[ match(common.mods, mods.historical)]

# useful
# https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
# start with strsplit, pmatch, charmatch

# ---------------------------------------------------------------------------------------
# Load historical model data
# ---------------------------------------------------------------------------------------

# Historical years
yrs = 1985:2014
nyr = length(yrs)

#variable matrix (turn this into a function?
varmat = matrix(nrow = length(fn.historical.kept), ncol = nyr)

for(i in 1:length(fn.historical.kept)){
# Open and extract data from the netcdf file
  fn = fn.historical.kept[i]
  fnpath = paste0("/data/users/hadaw/cmip6/areaavg/tas/", fn)

  nc = nc_open(fnpath)

  v = ncvar_get(nc, 'tas')
  nc_close(nc)
  
  vtail = tail(v, nyr)
  varmat[i, ] = vtail
}
  
matplot(t(varmat), type = 'l')

# Amomalize model runs to the mean of each run
varmat.means = apply(varmat, 1, mean)

# varmat.anom is the variable timeseries matrix anomaly
varmat.anom = sweep(varmat, 1, STATS = varmat.means)

# --------------------------------------------------------------------------------------------
# Load Observations
# --------------------------------------------------------------------------------------------

# Load up the HadCRUT4 data
obsmat = read.table(file = 'https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.annual_ns_avg.txt')

obsyears = obsmat[,1]
obsmed   = obsmat[,2]

obs.ix = which(obsyears %in% yrs)

obs = obsmed[obs.ix]
obs.anom = obs - mean(obs)

# Plot the timeseries and model historical runs.
matplot(yrs, t(varmat.anom), type = 'l', lty = 'solid', col = 'grey')
lines(yrs, obs.anom, col = 'red', lty = 'solid')
lines(yrs, varmat.anom[25,] ,col = 'blue', lty = 'solid')


# Matrix of tas errors (differences between the anomalies of the obs and models)
tas.err = sweep(varmat.anom, 2, obs.anom)

matplot(yrs, t(tas.err), type = 'l', lty = 'solid', col = 'grey')

tas.rmse = sqrt(apply( (tas.err)^2, 1, mean))

# Distance from observations for each of the models.
dotchart(tas.rmse, labels = common.mods)

# -------------------------------------------------------------------
# Weighting section from Sanderson et al. (2017)
# -------------------------------------------------------------------

# Skill weighting

D_q = 0.2 # guess for now

# w_q = exp(  - (tas.rmse/D_q)^2)

get_w_q = function(rmse, D_q){ 
  out = exp(-(rmse/D_q)^2)
  out
}

w_q = get_w_q(rmse = tas.rmse, D_q = D_q)

# A quick look at how D_q affects the skill weighting.
D_q.vec = seq(from = 0.1, to = 0.9, by = 0.1)

for(i in 1:length(D_q.vec)){

  w_q_i = get_w_q(rmse = tas.rmse, D_q = D_q.vec[i])
  tas.rmse.sort = sort(tas.rmse, index.return = TRUE)
  w_q_i.sort = w_q_i[tas.rmse.sort$ix]
  
  points(tas.rmse.sort$x, w_q_i.sort, type = 'o')

  text(
       tail(tas.rmse.sort$x,1),
       tail(w_q_i.sort,1),
       labels = D_q.vec[i],
       pos = 4,
       col = 'red'
       )
}
legend('topright', legend = 'D_q', text.col = 'red', bty = 'n')


# Similarity weighting (including observations)



# Similarity weighting
# (Not including the observations at the moment)

# measure pairwise distance between models and the obs
d_ij = dist(varmat.anom, upper = TRUE, diag = TRUE)
d_ij.matrix = as.matrix(d_ij)

# D_u is the radius of similarity
D_u = 0.5

# S is a similarity score
S = exp( - (d_ij.matrix / D_u)^2)

library(fields)

dev.new(width = 8, height = 8)

mat1 = apply(d_ij.matrix, 2, rev)

# I think this is the right way round.
dev.new(width = 10, height = 10)
par(mar = c(8,8,3,3))
image.plot(1:26, 1:26,t(mat1), axes = FALSE, xlab = '', ylab = '')
axis(2, at = 26:1, labels = common.mods, las = 1)
axis(1, at = 1:26, labels = common.mods, las = 2)


# R_u is a model's repetition
R = rep(NA, nrow(S))

for(i in 1:nrow(S)){
  r = sum(S[i,][-i])
  R[i] = r
}

R_u = 1 + R


# w_u is the Similarity weighting
w_u = 1/R_u

w.raw = w_q * w_u

# calculate A, a constant that ensures A sums to one
Ainv = sum(w.raw)
A = 1/Ainv

# final weights
w = w.raw*A


# I think that the independence weighting is dominating.
# Check how rmse and weighting interact


dev.new(width = 5, height = 8)
plot(tas.rmse, w_q, xlab = 'RMSE', ylab = 'weight', pty = 'n',
     ylim = c(0,1),
     xlim = c(min(tas.rmse), (max(tas.rmse)*1.05))
     )
points(tas.rmse, w_u, col = 'red')
points(tas.rmse, w, col = 'blue')
points(tas.rmse, w.raw, col = 'green')

                                        # here's how the similarity score and rmse relate.
plot(tas.rmse, w_u)





# Apply the weighting by sampling from the distribution of models, weighted
# according to w

nmods = nrow(varmat)

# create a much bigger timeseries matrix,

# create index by which we sample the model matrix
n = 100000
ix.w = sample(1:nmods, n, replace = TRUE, prob = w) # weighted sample
ix = sample(1:nmods, n, replace = TRUE)             # unweighted sample


ts.sample.weighted = varmat.anom[ix.w, ]
ts.sample.unweighted = varmat.anom[ix, ]

mat.change = function(X){
  cols = ncol(X)
  out = X[, cols] - X[,1]
  out
}
warming.weighted = mat.change(ts.sample.weighted)
warming.unweighted = mat.change(ts.sample.unweighted)


dev.new()
par(mfrow = c(2,1))
hist(warming.unweighted, xlim = c(0,1.5) )
hist(warming.weighted, xlim = c(0,1.5 ))

matplot(t(samp.tsmat), type = 'l')



# So, what is the mean warming over the unweighted
# and weighted sample?

mean(warming.weighted)
mean(warming.unweighted)


# --------------------------------------------------------------------------------------------
# Extract future projection data
#
#
# --------------------------------------------------------------------------------------------

# start with those models that are r1i1p1f1, to keep everything the same.
#fnallvec = dir('/data/users/hadaw/cmip6/areaavg/tas/', pattern = 'ssp585_r1i1p1')

yrs = 2015:2099
nyr = length(yrs)
futuremat = matrix(nrow = length(fn.ssp585.kept), ncol = nyr)

for(i in 1:length(fn.ssp585.kept)){
# Open and extract data from the netcdf file
  fn = fn.ssp585.kept[i]
  fnpath = paste0("/data/users/hadaw/cmip6/areaavg/tas/", fn)

  nc = nc_open(fnpath)

  v = ncvar_get(nc, 'tas')
  print(length(v))
  nc_close(nc)
  
  vhead = head(v, nyr)
 
  futuremat[i, ] = vhead
}



matplot(t(futuremat), type = 'l')

# Amomalize to the mean

# varmat.anom is the variable timeseries matrix anomaly
futuremat.changemat = sweep(futuremat, 1, STATS = futuremat[,1])
future.diff = futuremat.changemat [,ncol(futuremat.changemat )]

future.warming = sample(future.diff , n, replace = TRUE)
future.warming.weighted = sample(future.diff, n, replace = TRUE, prob = w)


dev.new()
par(mfrow = c(2,1))
hist(future.warming, xlim = c(0,10))
hist(future.warming.weighted, xlim  = c(0,10))


# I don't think my resampling technique works.

# try a very large number of samples for the mean 

mfw  = rep(NA,n)
mfww =  rep(NA,n)

for(i in 1:n){

  mfw[i] = mean(sample(future.diff, 26, replace = TRUE))
  mfww[i] = mean(sample(future.diff, 26, replace = TRUE, prob = w))
}
  
# same result




