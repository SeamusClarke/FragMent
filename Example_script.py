import FragMent
import numpy
import matplotlib.pyplot as plt
from scipy.stats import iqr
from astropy.io import fits
from matplotlib import gridspec

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
sep_limit = 0.06
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
nruns = 1000
lower_lim = 0.06
sep, w = FragMent.TwoPoint(pos,nruns,boundary,lower_lim)

### Calculate the approximate two point correlation function
op_crit = 1e-2
sep, wa, sig_w = FragMent.ApproxTwoPoint(pos,nruns,boundary,lower_lim,op_crit)

### Plot results
plt.figure(2)
plt.plot(sep,w,"b",label="Exact")
plt.plot(sep,wa,"k--",label="Approx")
plt.fill_between(sep,wa-sig_w,wa+sig_w,alpha=0.5)
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

### Place boundaries on the five parameters of the two-tier fragmentation model (the first 4 entries are the boundaries for the single-tier model)
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
n=15
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





### An example of straightening a filament and the accompanying cores

### Load the example filament map and spine
Map = fits.getdata("./Curved_filament.fits")
spine = numpy.loadtxt("./Curved_filament_spine.dat")

### Place some `cores' close to the spine
pos = numpy.zeros((4,2))

pos[0,0] = spine[5,0]
pos[0,1] = spine[5,1] + 1

pos[1,0] = spine[17,0] + 3
pos[1,1] = spine[17,1] + 1

pos[2,0] = spine[60,0]
pos[2,1] = spine[60,1]

pos[3,0] = spine[100,0] + 2
pos[3,1] = spine[100,1] - 3

### Map the cores to the spine
cpos = FragMent.Map_cores(spine,pos,10)

### Straighten the filament using the interpolation and weighted average methods
s,l,r = FragMent.Straighten_filament_interp(spine,Map.T,90,30,10)
s2,l,r = FragMent.Straighten_filament_weight(spine,Map.T,90,30,10,0.5)

### Get limits for the plots
miny = numpy.amin(l)
maxy = numpy.amax(l)

minx = numpy.amin(r)
maxx = numpy.amax(r)

maxc = numpy.amax(Map)
minc = numpy.amin(Map)

### Plot the results
fig = plt.figure(4,figsize=(10,8))
gs = gridspec.GridSpec(2,4)
ax1 = fig.add_subplot(gs[:,:2])

i = ax1.imshow(Map,origin=0,interpolation="gaussian", vmax = maxc, vmin = minc , extent=[0,200,0,200])
ax1.plot(spine[:,0],spine[:,1],"k",linewidth="2")
ax1.plot(pos[:,0],pos[:,1],"ro")
plt.setp(ax1,yticks=[40,80,120,160],yticklabels=[40,80,120,160],xticks=[40,80,120,160],xticklabels=[40,80,120,160])
plt.ylabel("y")
plt.xlabel("x")

ax2 = fig.add_subplot(gs[:,2])
ax2.imshow(s.T,origin="lower",interpolation="gaussian", vmax = maxc, vmin = minc , extent=[minx,maxx,miny,maxy])
ax2.plot(cpos[:,0],cpos[:,1],"ro")
plt.setp(ax2,yticks=[],yticklabels=[],xticks=[-15,0,15],xticklabels=[-15,0,15])
plt.xlabel("Radius")

ax3 = fig.add_subplot(gs[:,3])
ax3.imshow(s2.T,origin="lower",interpolation="gaussian", vmax = maxc, vmin = minc , extent=[minx,maxx,miny,maxy])
ax3.plot(cpos[:,0],cpos[:,1],"ro")
plt.setp(ax3,yticks=[20,40,60,80,100,120],yticklabels=[20,40,60,80,100,120],xticks=[-15,0,15],xticklabels=[-15,0,15])
plt.ylabel("Longitudinal axis")
plt.xlabel("Radius")
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")

cbaxes = fig.add_axes([0.1, 0.1, 0.8,0.05])
cb = plt.colorbar(i,label="Intensity", cax = cbaxes,orientation="horizontal") 

fig.subplots_adjust(left=0.11, bottom=0.18, right=0.89, top=0.93, wspace=0.25, hspace=0.1)

plt.savefig("Straight_filament.png")

