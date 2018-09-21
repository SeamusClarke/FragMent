import FragMent
import numpy
import matplotlib.pyplot as plt
from scipy.stats import iqr

### Load the positions of the cores
pos = numpy.loadtxt("core_positions.dat")

### Load the spine of the filaments
spine = numpy.loadtxt("spine.dat")

### Define a boundary for placing random cores
xmin = 1.994
xmax = 1.996
ymin = 0.0
ymax = 4.0
boundary = numpy.array([xmin,xmax,ymin,ymax])



### Perform a nearest neighbour test on the data
nn_seps = FragMent.NNS(pos)

### Statistical significance of the results
nruns = 1000
sep_limit = 0.02
p_median, p_mean = FragMent.Stat_Sig(nn_seps, "NNS", boundary, nruns, sep_limit)
ks_stat, p_ks = FragMent.KS_test(nn_seps, "NNS", boundary, nruns, sep_limit)
ad_stat, crit_vals, p_ad = FragMent.AD_test(nn_seps, "NNS", boundary, nruns, sep_limit)

### Report results
print " "
print " "
print "######## Nearest neighbour results ########"
print "The median and interquartile range of the distribution are: ", numpy.median(nn_seps), iqr(nn_seps)
print "The p-value using the median-interquartile range NHT is   : ", p_median
print "The mean and standard deviation of the distribution are   : ", numpy.mean(nn_seps), numpy.std(nn_seps)
print "The p-value using the mean-standard deviation NHT is      : ", p_mean
print "The p-values from the K-S and A-D test are                : ", p_ks, p_ad





### Perform a minimum spanning test on the data
mst_seps, mst = FragMent.MST(pos)

### Statistical significance of the results
p_median, p_mean = FragMent.Stat_Sig(mst_seps, "MST", boundary, nruns, sep_limit)
ks_stat, p_ks = FragMent.KS_test(mst_seps, "MST", boundary, nruns, sep_limit)
ad_stat, crit_vals, p_ad = FragMent.AD_test(mst_seps, "MST", boundary, nruns, sep_limit)

### Report results
print " "
print " "
print "######## Minimum spanning tree results ########"
print "The median and interquartile range of the distribution are: ", numpy.median(mst_seps), iqr(mst_seps)
print "The p-value using the median-interquartile range NHT is   : ", p_median
print "The mean and standard deviation of the distribution are   : ", numpy.mean(mst_seps), numpy.std(mst_seps)
print "The p-value using the mean-standard deviation NHT are     : ", p_mean
print "The p-values from the K-S and A-D test are                : ", p_ks, p_ad





### Perform a Nth nearest neighbour test on the data
NNN_sep = FragMent.NNNS(pos)

### Reduce the separation 2d-array into an average and a width
mean_seps = numpy.mean(NNN_sep, axis=0)
std_seps = numpy.std(NNN_sep, axis=0)
median_seps = numpy.median(NNN_sep, axis=0)
iqr_seps = iqr(NNN_sep, axis=0)

mean_seps = mean_seps[1:]
std_seps = std_seps[1:]
median_seps = median_seps[1:]
iqr_seps = iqr_seps[1:]

### Plot results
plt.figure(1)
plt.plot(mean_seps,std_seps,"b",label="Mean")
plt.plot(median_seps,iqr_seps,"r",label="Median")
plt.legend(loc="best")
plt.xlabel("Average separation")
plt.ylabel("Average width")
plt.savefig("NNN_results.png")

print " "
print " "
print "######## Nth nearest neighbour test ########"
print "The results are plotted in the file NNN_results.png"





### Determine the two-point correlation function

### Calculate the exact two point correlation function
nruns = 10000
lower_lim = 0.06
sep, w = FragMent.TwoPoint(pos,nruns,boundary,lower_lim)

### Calculate the approximate two point correlation function
op_crit = 1e-2
sep, wa = FragMent.ApproxTwoPoint(pos,nruns,boundary,lower_lim,op_crit)

### Plot results
plt.figure(2)
plt.plot(sep,w,"b",label="Exact")
plt.plot(sep,wa,"k--",label="Approx")
plt.xlabel("Separation")
plt.ylabel("Two point correlation function")
plt.legend(loc="best")
plt.savefig("Two_point_results.png")

print " "
print " "
print "######## Two point correlation function test ########"
print "The results are plotted in the file Two_point_results.png"


### Determine the power spectrum of the spine data
k,power = FragMent.FT_Spine(spine)

### Remove the k=0 mode
k = k[1:]
power = power[1:]

### Plot results
plt.figure(3)
plt.plot(k,power/power[0],"b")
plt.xlabel("Wavenumber")
plt.ylabel("Normalised power spectrum")
plt.savefig("FT_results.png")

print " "
print " "
print "######## Power spectrum test ########"
print "The results are plotted in the file FT_results.png"





### Frequentist model selection

### Place boundaries on the ten parameters of the two-tier fragmentation model (the first 4 are the boundaries for the single-tier model)
bound = numpy.zeros(10)

### Boundaries for the mean of the first Gaussian
bound[0] = 0.0
bound[1] = 1.0
### Boundaries for the standard deviation of the first Gaussian
bound[2] = 1e-99
bound[3] = 0.5
### Boundaries for the mean of the second Gaussian
bound[4] = 0.0
bound[5] = 1.0
### Boundaries for the standard deviation of the second Gaussian
bound[6] = 1e-99
bound[7] = 0.5
### Boundaries for the ratio of the amplitude of the two Gaussians
bound[8] = 0.0
bound[9] = 10.0

### Setup the MCMC parameter array
mcmc_param = numpy.zeros(3,dtype=int)
mcmc_param[0] = 100     ### Number of walkers
mcmc_param[1] = 200      ### Number of burn steps for each walker
mcmc_param[2] = 1000     ### Number of actual steps for each walker

### Run the model selection
model_one, theta_one, model_two, theta_two = FragMent.Freq_model_selection(mst_seps,bound,mcmc_param)

### Unpack outputs
aic_one = model_one[0]
weight_one = model_one[1]
aic_two = model_two[0]
weight_two = model_two[1]

### Report results
print " "
print " "
print "######## Frequentist model selection results ########"

if(aic_one < aic_two):
	print "Single-tier fragmentation is preferred"
	print "The Akaike weight for the model is: ", weight_one
	print "This gives an evidence ratio of   : ", weight_one/weight_two
	print "The fit parameters are            : ", theta_one

if(aic_two< aic_one):
	print "Two-tier fragmentation is preferred"
	print "The Akaike weight for the model is: ", weight_two
	print "This gives an evidence ratio of   : ", weight_two/weight_one
	print "The fit parameters are            : ", theta_two





### Bayesian model selection

### Declare the priors for model one and model two
### Both priors are flat, i.e. all model parameters have the same probability

def prior_one(params, bound):

	mu1, sig1 = params

	term1 = 1.0/(bound[1]-bound[0])
	term2 = 1.0/(bound[3]-bound[2])

	return term1 * term2

def prior_two(params, bound):

	mu1, sig1, mu2, sig2, ratio = params

	term1 = 1.0/(bound[1]-bound[0])
	term2 = 1.0/(bound[3]-bound[2])
	term3 = 1.0/(bound[5]-bound[4])
	term4 = 1.0/(bound[7]-bound[6])
	term5 = 1.0/(bound[9]-bound[8])

	return term1*term2*term3*term4*term5


### Calculate the posterior for each model, these can be plotting using corner
walkers_one, post_one, walkers_two, post_two = FragMent.Posterior_of_models(mst_seps,bound,mcmc_param,prior_one, prior_two)

### Calculate the evidence of each model
n=25
evi_one, evi_two = FragMent.Evidence_explicit(mst_seps,bound,prior_one, prior_two,n)

### Report results
print " "
print " "
print "######## Bayesian model selection results ########"

if(evi_one > evi_two):
	print "Single-tier fragmentation is preferred"
	print "The Bayes factor for in favour of model one is   : ", evi_one/evi_two
	print "The most likely parameters from the posterior are: ", walkers_one[numpy.argmax(post_one)]

if(evi_two > evi_one):
	print "Two-tier fragmentation is preferred"
	print "The Bayes factor for in favour of model two is   : ", evi_two/evi_one
	print "The most likely parameters from the posterior are: ", walkers_two[numpy.argmax(post_two)]

