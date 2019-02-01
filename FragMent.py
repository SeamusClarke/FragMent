### Import numpy and the 1D fft function
import numpy
from numpy.fft import rfft

### Import the scipy parts needed for minimum spanning trees
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

### Import scipy stats and commonly used functions
import scipy.stats
from scipy.stats import iqr
from scipy.stats import gaussian_kde
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp

### Library for python to C communication
import ctypes
from numpy.ctypeslib import ndpointer

### emcee for MCMCs
import emcee



########################################################################################
################### Function which are called from outside FragMent ####################
########################################################################################



### Statistical significance tools ###



### A function to perform an Anderson-Darling test on the results from either the nearest neighbour method or the minimum spanning tree
### Inputs:
### dist         (array)                   The distribution of separations which is being tested
### method       (string)                  The method used to produce the distribution of separations
### bounds        (array)                   The boundary box in which one may randomly place cores
### num          (int)                     The number of random realisations to produce the null distribution
### lower_lim    (float)                   The smallest separation between cores allowed in the random distributions, should be set to data resolution.
###
### Outputs:
### stat_dist    (float)                   The AD test statistic
### crit_val     (array)                   An array showing the critical values of the test statistic to achieve a significance levels 25%, 10%, 5%, 2.5%, 1%
### p_val        (float)                   The p-value corresponding to the test statistic 

def AD_test(dist,method,bounds,num,lower_lim):

	### Unpack the bounds array
	xmin = bounds[0]
	xmax = bounds[1]

	ymin = bounds[2]
	ymax = bounds[3]

	### Finding the number of cores to randomly place from the input distribution
	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	### Initialise a random seed
	numpy.random.seed(1)

	### Create counter and array for the separations from randomly placed cores 
	tot_num = 0
	tot_sep = numpy.array([],dtype=float)

	### Loop over each random realisation
	while(tot_num<num):

		### Place the cores randomly in the box
		x = (xmax - xmin)*numpy.random.random(num_cores) + xmin
		y = (ymax - ymin)*numpy.random.random(num_cores) + ymin
		pos = numpy.column_stack((x,y))

		### Use the same method to analyse the random cores as the data set
		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)

		### Check that no randomly placed cores are too close to each other
		if(len(seps)!=numpy.sum(seps>lower_lim)):
			continue

		### If the cores are fine add them to the array for the null distribution
		tot_num = tot_num + 1
		tot_sep = numpy.concatenate((tot_sep,seps))

	### Perform the AD test between the data set distribution and the null distribution
	stat_dist,crit_val,p_val = anderson_ksamp([dist,tot_sep])
	return stat_dist,crit_val,p_val











### A function to perform a Kolmogorov-Smirnov test on the results from either the nearest neighbour method or the minimum spanning tree
### Inputs:
### dist         (array)                   The distribution of separations which is being tested
### method       (string)                  The method used to produce the distribution of separations
### bounds        (array)                   The boundary box in which one may randomly place cores
### num          (int)                     The number of random realisations to produce the null distribution
### lower_lim    (float)                   The smallest separation between cores allowed in the random distributions, should be set to data resolution.
###
### Outputs:
### stat_dist    (float)                   The KS test statistic
### p_val        (float)                   The p-value corresponding to the test statistic

def KS_test(dist,method,bounds,num,lower_lim):

	### Unpack the bounds array
	xmin = bounds[0]
	xmax = bounds[1]

	ymin = bounds[2]
	ymax = bounds[3]

	### Finding the number of cores to randomly place from the input distribution
	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	### Initialise a random seed
	numpy.random.seed(1)

	### Create counter and array for the separations from randomly placed cores
	tot_num = 0
	tot_sep = numpy.array([],dtype=float)

	### Loop over each random realisation
	while(tot_num<num):

		### Place the cores randomly in the box
		x = (xmax - xmin)*numpy.random.random(num_cores) + xmin
		y = (ymax - ymin)*numpy.random.random(num_cores) + ymin
		pos = numpy.column_stack((x,y))

		### Use the same method to analyse the random cores as the data set
		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)

		### Check that no randomly placed cores are too close to each other
		if(len(seps)!=numpy.sum(seps>lower_lim)):
			continue

		### If the cores are fine add them to the array for the null distribution
		tot_num = tot_num + 1
		tot_sep = numpy.concatenate((tot_sep,seps))


	### Perform the KS test between the data set distribution and the null distribution
	stat_dist,p_val = ks_2samp(dist,tot_sep)
	return stat_dist,p_val










### A function to perform a null hypothesis test using the average and width of the distribution 
### Inputs:
### dist         (array)                   The distribution of separations which is being tested
### method       (string)                  The method used to produce the distribution of separations
### bound        (array)                   The boundary box in which one may randomly place cores
### num          (int)                     The number of random realisations to produce the null distribution
### lower_lim    (float)                   The smallest separation between cores allowed in the random distributions, should be set to data resolution.
###
### Outputs:
### p_med        (float)                   The p-value from the null hypothesis test when using the median and interquartile range
### p_mea        (float)                   The p-value from the null hypothesis test when using the mean and standard deviation

def Stat_Sig(dist, method, bounds, num, lower_lim):

	### Unpack the bounds array
	xmin = bounds[0]
	xmax = bounds[1]

	ymin = bounds[2]
	ymax = bounds[3]

	### Create empty lists for the medians, means, standard deviation and interquartile ranges of the randomly placed core distributions
	med_sep = []
	mea_sep = []
	std_sep = []
	iqr_sep = []

	### Finding the number of cores to randomly place from the input distribution	
	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	### Initialise a random seed
	numpy.random.seed(1)

	### Create counter for the number of random realisations
	tot_num = 0

	### Loop over each random realisation
	while(tot_num<num):

		### Place the cores randomly in the box
		x = (xmax - xmin)*numpy.random.random(num_cores) + xmin
		y = (ymax - ymin)*numpy.random.random(num_cores) + ymin
		pos = numpy.column_stack((x,y))

		### Use the same method to analyse the random cores as the data set
		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)

		### Check that no randomly placed cores are too close to each other
		if(len(seps)!=numpy.sum(seps>lower_lim)):
			continue

		### Add the average and width measurements for this random realisation to the storage list
		med_sep.append(numpy.median(seps))
		iqr_sep.append(iqr(seps))

		mea_sep.append(numpy.mean(seps))
		std_sep.append(numpy.std(seps))

		### Update counter
		tot_num = tot_num + 1

	### Transform lists to arrays and then stack them
	med_sep = numpy.array(med_sep,dtype=float)
	mea_sep = numpy.array(mea_sep,dtype=float)

	iqr_sep = numpy.array(iqr_sep,dtype=float)
	std_sep = numpy.array(std_sep,dtype=float)

	med_iqr = numpy.column_stack((med_sep,iqr_sep))
	mea_std = numpy.column_stack((mea_sep,std_sep))

	### Find the average and width of the input distribution
	med_dist = numpy.median(dist)
	iqr_dist = iqr(dist)

	mea_dist = numpy.mean(dist)
	std_dist = numpy.std(dist)



	### Determine the maximum of the median and interquartile range to produce limits for the KDE
	maxx = numpy.amax(med_sep)
	if(med_dist>maxx):
		maxx = med_dist

	maxy = numpy.amax(iqr_sep)
	if(iqr_dist>maxy):
		maxy = iqr_dist
	
	### Produce the KDE for the median and interquartile range and produce a 100 by 100 square of evaluation points 
	kde_med = gaussian_kde(med_iqr.T,bw_method="scott")
	x1,y1 = numpy.mgrid[0:1.1*maxx:100j,0:1.1*maxy:100j]
	po = numpy.vstack([x1.ravel(),y1.ravel()])

	### Evaluate the KDE at these positions
	k_med = numpy.reshape(kde_med(po).T,x1.shape)
	k_med = k_med.T

	### Determine the position of the data-set's median and interquartile range in this median-interquartile range space 
	maxx = numpy.amax(x1)
	minx = numpy.amin(x1)
	dx = (maxx - minx)/ 100.
	
	maxy = numpy.amax(y1)
	miny = numpy.amin(y1)
	dy = (maxy - miny)/ 100.

	ix = int((med_dist-minx)/dx)
	iy = int((iqr_dist-miny)/dy)

	### Determine the probability of the data-sets median and interquartile range given the null hypothesis and the positions which are more extreme
	val = k_med[iy,ix]
	less = k_med[k_med<=val]

	### The p-value is calculated by summing all outcomes as, or more, extreme than the data-set and normalising over all possible outcomes.
	p_med = numpy.sum(less) / numpy.sum(k_med)



	### Determine the maximum of the mean and standard deviation to produce limits for the KDE
	maxx = numpy.amax(mea_sep)
	if(mea_dist>maxx):
		maxx = mea_dist

	maxy = numpy.amax(std_sep)
	if(std_dist>maxy):
		maxy = std_dist

	### Produce the KDE for the mean and standard deviation and produce a 100 by 100 square of evaluation points 
	kde_mea = gaussian_kde(mea_std.T,bw_method="scott")
	x2,y2 = numpy.mgrid[0:1.1*maxx:100j,0:1.1*maxy:100j]
	po2 = numpy.vstack([x2.ravel(),y2.ravel()])

	### Evaluate the KDE at these positions
	k_mea = numpy.reshape(kde_mea(po2).T,x2.shape)
	k_mea = k_mea.T

	### Determine the position of the data-set's mean and standard deviation in this mean-standard deviation space 
	maxx = numpy.amax(x2)
	minx = numpy.amin(x2)
	dx = (maxx - minx)/ 100.
	
	maxy = numpy.amax(y2)
	miny = numpy.amin(y2)
	dy = (maxy - miny)/ 100.

	ix = int((mea_dist-minx)/dx)
	iy = int((std_dist-miny)/dy)

	### Determine the probability of the data-sets mean and standard deviation given the null hypothesis and the positions which are more extreme
	val = k_mea[iy,ix]
	less = k_mea[k_mea<=val]

	### The p-value is calculated by summing all outcomes as, or more, extreme than the data-set and normalising over all possible outcomes.
	p_mea = numpy.sum(less) / numpy.sum(k_mea)



	### Return the p-values for the median-interquartile and mean-standard deviation measures.
	return p_med, p_mea










### Fragmentation analysis techniques ###



### A function to produce the nearest neighbour separation distribution for a given set of core positions 
### Inputs:
### pos        (2d-array)                The 2d-array of size (n,2) for n cores 
###
### Outputs:
### seps       (array)                   The distribution of separations to each core's nearest neighbour

def NNS(pos):

	### Determine the number of cores
	n_pos = len(pos[:,0])

	### Create an empty array for the separations
	seps = numpy.zeros(n_pos)

	### Loop over all cores
	for ii in range(0,n_pos):

		### nearest-neighbour separation is initialised with a large value
		min_sep = 1e99

		### Loop over all cores to find the separations
		for jj in range(0,n_pos):

			### Skip if determining the separation between a core and itself
			if(ii==jj):
				continue

			### Determine the separation between the two cores
			sep = numpy.sqrt( (pos[ii,0]-pos[jj,0])**2 + (pos[jj,1]-pos[ii,1])**2 )

			### If the separation is the smaller than the minimum separation so far then replace the minimum separation
			if(sep<min_sep):
				min_sep = sep

		### Once all separations for this core have been considered store the minimum separation in the nearest neighbour separation array
		seps[ii] = min_sep

	### Return the distribution of separations
	return seps










### A function to produce the Nth nearest neighbour separation distribution for a given set of core positions 
### Inputs:
### pos        (2d-array)               The 2d-array of size (n,2) for n cores 
###
### Outputs:
### sep       (2d-array)                The 2d-array of size (n,n) which gives the ordered separation distribution between every pair of cores

def NNNS(pos):

	### Determine the number of cores
	n_pos = len(pos[:,0])

	### Create an empty array for the separations
	sep = numpy.zeros((n_pos,n_pos))

	### Loop over all cores
	for ii in range(0,n_pos):

		### Loop over all cores
		for jj in range(0,n_pos):

			### Calculate the separation between the two cores
			sep[ii,jj] = numpy.sqrt( (pos[ii,0]-pos[jj,0])**2 + (pos[jj,1]-pos[ii,1])**2 )

		### Sort the separation distribution array so that sep[:,N] will show the Nth nearest neighbour separations for all n cores 
		sep[ii,:] = numpy.sort(sep[ii,:])

	### Return the separation distribution array
	return sep










### A function to produce minimum spanning tree for a set of core positions and the distribution of resulting edge lengths 
### Inputs:
### pos        (2d-array)               The 2d-array of size (n,2) for n cores 
###
### Outputs:
### dis        (array)                  The distribution of edge lengths from the minimum spanning tree
### mst        (array)                  The minimum spanning tree as a matrix

def MST(pos):

	### Determine the number of cores
	n_pos = len(pos[:,0])
	
	### Array to store the complete graph
	DD = numpy.zeros((n_pos,n_pos))

	### Loop over all cores to determine the complete graph
	for jj in range(0,n_pos):
		for kk in range(jj,n_pos):

			### Store the separation between every core pair in the complete graph taking advantage of the fact it is symmetric
			DD[jj,kk] = numpy.sqrt( (pos[jj,0]-pos[kk,0])**2 + (pos[jj,1]-pos[kk,1])**2 )
			DD[kk,jj] = DD[jj,kk]

	### Turn the complete graph into a sparse matrix 
	DD2 = csr_matrix(DD)

	### Construct the minimum spanning tree from the sparse matrix complete graph and transform it into an array
	t = minimum_spanning_tree(DD2)
	mst = t.toarray().astype(float)

	### Flatten the 2d mst array into a 1d array and remove all the zeros
	dis = mst.flatten()
	dis = dis[dis>0]

	### Return the distribution of edge lengths and the minimum spanning tree array for plotting
	return numpy.array(dis,dtype=float), mst










### A function to produce the (exact) two-point correlation function for a given set of core positions 
### Inputs:
### pos        (2d-array)               The 2d-array of size (n,2) for n cores 
### nruns      (int)                    The number of random realisations to produce the DR and RR arrays
### bounds     (array)                  The boundary box in which one may randomly place cores
### lower_lim  (float)                  The smallest value of the bandwidth allowed
###
### Outputs:
### sep       (array)                   The array of separations at which the two-point correlation function has been evaluated
### w         (array)                   The two-point correlation function 

def TwoPoint(pos,nruns,bounds,lower_lim):

	### Determine the number of cores
	n_pos = len(pos[:,0])

	### Unpack the positions into two separate x and y arrays and ensure they are doubles
	pos_x = pos[:,0]
	pos_y = pos[:,1]

	pos_x = numpy.array(pos_x,dtype=numpy.double)
	pos_y = numpy.array(pos_y,dtype=numpy.double)

	### Unpack the bounds array
	xmin = bounds[0]
	xmax = bounds[1]

	ymin = bounds[2]
	ymax = bounds[3]

	### Set the number of evaluation points and produce the array for these evaluation points
	x_size=1000
	max_dist = numpy.sqrt((xmax-xmin)**2+(ymax-ymin)**2)
	sep = numpy.linspace(0,max_dist,x_size)

	### Create empty arrays for the DD, DR and RR KDEs 
	DD = numpy.zeros(x_size, dtype=numpy.double)
	DR = numpy.zeros(x_size, dtype=numpy.double)
	RR = numpy.zeros(x_size, dtype=numpy.double)

	### Communication between python and the underlying C code
	lib = ctypes.cdll.LoadLibrary("./FragMent_C.so")
	twopoint = lib.TwoPoint
	twopoint.restype = None
	lower_lim = numpy.double(lower_lim)
	twopoint.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_double, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

	### Call the C two-point correlation function which returns KDEs for the DD, DR and RR arrays
	twopoint(pos_x,pos_y,n_pos,nruns,bounds,x_size,lower_lim,DD,DR,RR)

	### Construct the two-point correlation function from the DD, DR and RR arrays
	w = (DD - 2*DR + RR)/RR	

	### Return the two-point correlation function and the array of separations at which it has been evaluated
	return sep, w










### A function to produce the (approximate) two-point correlation function for a given set of core positions 
### Inputs:
### pos        (2d-array)               The 2d-array of size (n,2) for n cores 
### nruns      (int)                    The number of random realisations to produce the DR and RR arrays
### bounds     (array)                  The boundary box in which one may randomly place cores
### lower_lim  (float)                  The smallest value of the bandwidth allowed
### op_crit    (op_crit)                The opening criterion constant to control the error of the approximate KDE
###
### Outputs:
### sep       (array)                   The array of separations at which the two-point correlation function has been evaluated
### w         (array)                   The two-point correlation function 
### sigma_w   (array)                   The error on the two-point correlation function (Taken from Grazian et al. 2006, A&A, 453, 507-515)

def ApproxTwoPoint(pos,nruns,bounds,lower_lim,error):

	### Determine the number of cores
	n_pos = len(pos[:,0])

	### Unpack the positions into two separate x and y arrays and ensure they are doubles
	pos_x = pos[:,0]
	pos_y = pos[:,1]

	pos_x = numpy.array(pos_x,dtype=numpy.double)
	pos_y = numpy.array(pos_y,dtype=numpy.double)

	### Unpack the bounds array
	xmin = bounds[0]
	xmax = bounds[1]

	ymin = bounds[2]
	ymax = bounds[3]

	### Set the number of evaluation points and produce the array for these evaluation points
	x_size=1000
	max_dist = numpy.sqrt((xmax-xmin)**2+(ymax-ymin)**2)
	sep = numpy.linspace(0,max_dist,x_size)

	### Create empty arrays for the DD, DR and RR KDEs 
	DD = numpy.zeros(x_size, dtype=numpy.double)
	DR = numpy.zeros(x_size, dtype=numpy.double)
	RR = numpy.zeros(x_size, dtype=numpy.double)

	### Communication between python and the underlying C code
	lib = ctypes.cdll.LoadLibrary("./FragMent_C.so")
	twopoint = lib.ApproxTwoPoint
	twopoint.restype = None
	lower_lim = numpy.double(lower_lim)
	twopoint.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_double, ctypes.c_double, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

	### Call the C two-point correlation function which returns KDEs for the DD, DR and RR arrays
	twopoint(pos_x,pos_y,n_pos,nruns,bounds,x_size,lower_lim,error,DD,DR,RR)
	
	### Construct the two-point correlation function from the DD, DR and RR arrays
	w = (DD - 2*DR + RR)/RR	

	### Determine the error on the two-point correlation function
	sigma_w = numpy.sqrt((1+w)/(n_pos*(n_pos-1)*DD+1e-10))

	### Return the two-point correlation function and the array of separations at which it has been evaluated
	return sep, w, sigma_w










### A function to produce a power spectrum of the Fourier transform of the filament's spine profile
### Inputs:
### spine      (array)                  The spine of the filament (i.e. column density / integrated intensity along the spine)
###
### Outputs:
### k          (array)                  The wavenumber array 
### Power      (array)                  The power spectrum 

def FT_Spine(spine):

	### Produce a Fourier transform of the spine
	fspine = rfft(spine)

	### Produce the wavenumber array, running from 0 to n-1 for a spine array of length n
	k = numpy.linspace(0,len(fspine)-1,len(fspine))

	### Take the square of the absolute of the Fourier transform to determine the power spectrum
	Power = numpy.abs(fspine)**2

	### Return the wavenumber and power spectrum	
	return k,Power

	








### Model selection techniques ###



### A function to determine the Akaike information criterion and other information to select between the single-tiered (model one) and two-tiered (model two) model
### Inputs:
### spacings      (array)               The distribution of edge lengths resulting from applying the minimum spanning tree technique on the data
### bound         (array)               The boundaries for the model parameters
### mcmc_params   (array)               An array containing the number of walkers, number of burn steps and number of "real" steps for the MCMC
###
### Outputs:
### model_one     (array)               An array containing the AICc, Akaike weight and log maximum likelihood for model one
### theta_one     (array)               An array containing the parameters for model one which result in the maximum likelihood
### model_two     (array)               An array containing the AICc, Akaike weight and log maximum likelihood for model two
### theta_two     (array)               An array containing the parameters for model two which result in the maximum likelihood

def Freq_model_selection(spacings,bound,mcmc_params):

	### Unpack the MCMC parameter array
	NWALKERS = mcmc_params[0]
	NBURN = mcmc_params[1]
	NSTEPS = mcmc_params[2]

	### Initialise the walkers for model one with random values within the parameter bounds
	NDIMS = 2
	theta0 = numpy.zeros((NWALKERS,NDIMS))
	theta0[:,0] = numpy.random.uniform(bound[0],bound[1],NWALKERS)
	theta0[:,1] = numpy.random.uniform(bound[2],bound[3],NWALKERS)	

	### Using emcee's to initialise the walkers
	sampler = emcee.EnsembleSampler(NWALKERS, NDIMS, lnlike_one_gauss, args=[spacings, bound[:4]])

	### Run burn steps
	theta, prob, state = sampler.run_mcmc(theta0,NBURN)

	### Run actual steps starting from the last of the burn values
	theta2, prob2, state2 = sampler.run_mcmc(theta, NSTEPS)

	### Unpack all of the walkers from the MCMC 
	walks_one = sampler.flatchain

	### Determine the log-likelihood for each value along the walkers 
	lnlike_one = numpy.zeros(len(walks_one[:,0]))
	for ii in range(0,len(lnlike_one)):		
		lnlike_one[ii] = lnlike_one_gauss(walks_one[ii,:], spacings, bound)
	
	### Determine the set of parameters which has the largest log-likelihood
	theta_one = walks_one[numpy.argmax(lnlike_one),:]

	### Store the maximum log-likelihood
	prob_one = numpy.amax(lnlike_one)

	### Work out the AICc for model one
	aic_one = 4 - 2*prob_one + 12.0/(len(spacings)-3 + 1e-99)

	

	### Initialise the walkers for model two with random values within the parameter bounds
	NDIMS = 5
	theta0 = numpy.zeros((NWALKERS,NDIMS))
	theta0[:,0] = numpy.random.uniform(bound[0],bound[1],NWALKERS)
	theta0[:,1] = numpy.random.uniform(bound[2],bound[3],NWALKERS)
	theta0[:,2] = numpy.random.uniform(bound[4],bound[5],NWALKERS)
	theta0[:,3] = numpy.random.uniform(bound[6],bound[7],NWALKERS)
	theta0[:,4] = numpy.random.uniform(bound[8],bound[9],NWALKERS)
	
	### Using emcee to initialise the walkers
	sampler = emcee.EnsembleSampler(NWALKERS, NDIMS, lnlike_two_gauss, args=[spacings, bound])

	### Run burn steps
	theta, prob, state = sampler.run_mcmc(theta0,NBURN)

	### Run actual steps starting from the last of the burn values
	theta2, prob2, state2 = sampler.run_mcmc(theta, NSTEPS)
	
	### Unpack all of the walkers from the MCMC 
	walks_two = sampler.flatchain

	### Determine the log-likelihood for each value along the walkers 
	lnlike_two = numpy.zeros(len(walks_two[:,0]))
	for ii in range(0,len(lnlike_one)):
		lnlike_two[ii] = lnlike_two_gauss(walks_two[ii,:], spacings, bound)

	### Determine the set of parameters which has the largest log-likelihood
	theta_two = walks_two[numpy.argmax(lnlike_two),:]

	### Store the maximum log-likelihood
	prob_two = numpy.amax(lnlike_two)

	### Work out the AICc for model two
	aic_two = 10 - 2*prob_two + 60.0/(len(spacings)-6 + 1e-99)



	### Determine the model with the minimum AICc and thus compute the Akaike weights
	min_aic = numpy.amin([aic_one,aic_two])
	da_one = abs(aic_one - min_aic)
	da_two = abs(aic_two - min_aic)

	weight_one = numpy.exp(-0.5*da_one) / (numpy.exp(-0.5*da_one) + numpy.exp(-0.5*da_two))
	weight_two = numpy.exp(-0.5*da_two) / (numpy.exp(-0.5*da_one) + numpy.exp(-0.5*da_two))

	### Package the information into output arrays
	model_one = numpy.array([aic_one,weight_one],dtype=float)
	model_two = numpy.array([aic_two,weight_two],dtype=float)

	### Return outputs
	return model_one, theta_one, model_two, theta_two










### A function to return un-normalised posteriors for the single-tiered (model one) and two-tiered (model two) models
### Inputs:
### spacings      (array)               The distribution of edge lengths resulting from applying the minimum spanning tree technique on the data
### bound         (array)               The boundaries for the model parameters
### mcmc_params   (array)               An array containing the number of walkers, number of burn steps and number of "real" steps for the MCMC
### prior_one     (function)            The prior for the model one parameters
### prior_two     (function)            The prior for the model two parameters
###
### Outputs:
### walks_one     (array)               An array containing the parameters at all walker positions for model one
### post_one      (array)               An array containing the posterior value at those walker positions for model one
### walks_two     (array)               An array containing the parameters at all walker positions for model two
### post_two      (array)               An array containing the posterior value at those walker positions for model two

def Posterior_of_models(spacings,bound,mcmc_params,prior_one, prior_two):

	### Unpack the MCMC parameter array
	NWALKERS = mcmc_params[0]
	NBURN = mcmc_params[1]
	NSTEPS = mcmc_params[2]

	### Initialise the walkers for model one with random values within the parameter bounds
	NDIMS = 2
	theta0 = numpy.zeros((NWALKERS,NDIMS))
	theta0[:,0] = numpy.random.uniform(bound[0],bound[1],NWALKERS)
	theta0[:,1] = numpy.random.uniform(bound[2],bound[3],NWALKERS)	

	### Using emcee to initialise the walkers
	sampler = emcee.EnsembleSampler(NWALKERS, NDIMS, post_one_gauss, args=[spacings, bound[:4], prior_one])

	### Run burn steps
	theta, prob, state = sampler.run_mcmc(theta0,NBURN)

	### Run actual steps starting from the last of the burn values
	theta2, prob2, state2 = sampler.run_mcmc(theta, NSTEPS)

	### Unpack all of the walkers from the MCMC 
	walks_one = sampler.flatchain

	### Determine the posterior for each value along the walkers 
	post_one = numpy.zeros(len(walks_one[:,0]))
	for ii in range(0,len(post_one)):
		post_one[ii] = post_one_gauss(walks_one[ii,:], spacings, bound[:4], prior_one)
	


	### Initialise the walkers for model two with random values within the parameter bounds
	NDIMS = 5
	theta0 = numpy.zeros((NWALKERS,NDIMS))
	theta0[:,0] = numpy.random.uniform(bound[0],bound[1],NWALKERS)
	theta0[:,1] = numpy.random.uniform(bound[2],bound[3],NWALKERS)
	theta0[:,2] = numpy.random.uniform(bound[4],bound[5],NWALKERS)
	theta0[:,3] = numpy.random.uniform(bound[6],bound[7],NWALKERS)
	theta0[:,4] = numpy.random.uniform(bound[8],bound[9],NWALKERS)

	### Using emcee to initialise the walkers	
	sampler = emcee.EnsembleSampler(NWALKERS, NDIMS, post_two_gauss, args=[spacings, bound, prior_two])

	### Run burn steps
	theta, prob, state = sampler.run_mcmc(theta0,NBURN)

	### Run actual steps starting from the last of the burn values
	theta2, prob2, state2 = sampler.run_mcmc(theta, NSTEPS)
	
	### Unpack all of the walkers from the MCMC 
	walks_two = sampler.flatchain

	### Determine the posterior for each value along the walkers 
	post_two = numpy.zeros(len(walks_two[:,0]))
	for ii in range(0,len(post_two)):	
		post_two[ii] = post_two_gauss(walks_two[ii,:], spacings, bound, prior_two)


	
	### Return the outputs
	return walks_one, post_one, walks_two, post_two
	
	








### A function to return the Bayesian evidence for the single-tiered (model one) and two-tiered (model two) models
### Inputs:
### spacings      (array)               The distribution of edge lengths resulting from applying the minimum spanning tree technique on the data
### bound         (array)               The boundaries for the model parameters
### prior_one     (function)            The prior for the model one parameters
### prior_two     (function)            The prior for the model two parameters
### n             (int)                 The number of evaluation points along a single axis to determine the evidence
###
### Outputs:
### evi_one       (float)               The evidence for model one
### evi_two       (float)               The evidence for model two

def Evidence_explicit(spacings,bound,prior_one, prior_two,n):
	
	### Determine the spacing between the evaluation points for model one
	dm = (bound[1]-bound[0])/float(n)
	ds = (bound[3]-bound[2])/float(n)

	### Create array to store the posterior
	post_one = numpy.zeros((n,n))

	### Loop over all evaluation points to determine the posterior
	for ii in range(0,n):
		for jj in range(0,n):

			### Determine the model parameters at this evaluation point
			mean = bound[0] + dm*ii
			sig = bound[2] + ds*jj
			params = [mean,sig]

			### Determine and store the posterior of model one at each evaluation point
			post_one[ii,jj] = post_one_gauss(params, spacings, bound, prior_one)

	### Produce a weight array so that the integral of the un-normalised posterior can be determined using the trapezium rule
	w = numpy.ones_like(post_one)
	w[0,:] = 0.5
	w[-1,:] = 0.5
	w[:,0] = 0.5
	w[:,-1] = 0.5

	### Determine the evidence of model one
	evi_one = numpy.sum(w*post_one)*dm*ds



	### Determine the spacing between the evaluation points for model two
	dm2 = (bound[5]-bound[4])/float(n)
	ds2 = (bound[7]-bound[6])/float(n)
	dr  = (bound[9]-bound[8])/float(n)

	### Create array to store the posterior
	post_two = numpy.zeros((n,n,n,n,n))

	### Loop over all evaluation points to determine the posterior
	for ii in range(0,n):
		for jj in range(0,n):
			for kk in range(0,n):
				for ll in range(0,n):
					for mm in range(0,n):

						### Determine the model parameters at this evaluation point
						mean1 = bound[0] + dm*ii
						sig1 = bound[2] + ds*jj
						mean2 = bound[4] + dm2*kk
						sig2 = bound[6] + ds2*ll
						r = bound[8] + dr*mm
						params = [mean1,sig1,mean2,sig2,r]

						### Determine and store the posterior of model two at each evaluation point
						post_two[ii,jj,kk,ll,mm] = post_two_gauss(params, spacings, bound, prior_two)

	### Produce a weight array so that the integral of the un-normalised posterior can be determined using the trapezium rule
	w = numpy.ones_like(post_two)
	w[0,:,:,:,:] = 0.5
	w[-1,:,:,:,:] = 0.5
	w[:,0,:,:,:] = 0.5
	w[:,-1,:,:,:] = 0.5
	w[:,:,0,:,:] = 0.5
	w[:,:,-1,:,:] = 0.5
	w[:,:,:,0,:] = 0.5
	w[:,:,:,-1,:] = 0.5
	w[:,:,:,:,0] = 0.5
	w[:,:,:,:,-1] = 0.5

	### Determine the evidence of model two
	evi_two = numpy.sum(w*post_two)*dm*ds*dm2*ds2*dr



	### Return the evidence of model one and model two
	return evi_one, evi_two

	








### Supplementary functions ###


### A function to straighten a filament using the weighted average method 
### Inputs:
### spine        (2d-array)             An ordered list of x,y co-ordinates of the spine points, shape n by 2 for n spine points
### Map          (2d-array)             The column density or integrated intensity map of the filament which you wish to straighten
### n_pix        (int)                  The number of evaluation points on either side of the radial profile
### max_dist     (float)                The distance to which the radial profile extends
### order        (int)                  The order of polynomial used to fit the spine points
### h_length     (float)                The smoothing length used for weighting       
###
### Str_fil      (2d-array)             The straightened filament map
### length       (array)                The length of the straightened filament
### radial       (array)                The radius at each evaluation point

def Straighten_filament_weight(spine, Map, n_pix, max_dist, order, h_length): 

	### Unpack the spine array and determine the number of spine points
	x = spine[:,0]
	y = spine[:,1]
	n_spine = len(x)

	### Fit a polynomial to the spine points, then find the gradient by differentiating the function 
	p = numpy.polyfit(x,y,order)
	f = numpy.poly1d(p)
	p2 = p[:-1]
	q = numpy.arange(1,order+1,1)
	grad = numpy.poly1d(p2*q[::-1])

	### Create an empty array for the straightened filament
	Str_fil = numpy.zeros((2*n_pix + 1,n_spine))

	### Pixels which are within 3 times the smoothing length are considered
	eval_dist = int(numpy.ceil(3*h_length)) + 1

	### Create empty length array
	length = numpy.zeros(n_spine)

	### Loop over the spine points
	for ii in range(0,n_spine):

		### Store co-ordinates of the spine point
		x0 = int(x[ii])
		y0 = int(y[ii])

		### If the spine point is not the first then we construct the length along the filament
		if(ii>0):
			length[ii] = length[ii-1] + numpy.sqrt( (x[ii]-x[ii-1])**2 + (y[ii]-y[ii-1])**2 )

		### Find the gradient at this spine point and using it as a normal define the radial axis perpendicular to the spine
		a = 1
		b = grad(x0)
		norm = numpy.array([a,b],dtype=numpy.double)
		plane_const = numpy.sum(norm*spine[ii,:]) 

		### Construct the evaluation points along the normal line
		m = -a/b
		xmax = x0 + numpy.sqrt(max_dist**2/(m*m + 1))
		xmin = x0 - numpy.sqrt(max_dist**2/(m*m + 1))
		x_norm = numpy.linspace(xmin,xmax,2*n_pix + 1)
		y_norm = m*x_norm - m*x0 + y0

		### Set up profile
		profile = numpy.zeros(2*n_pix + 1)
		profile_weight = numpy.zeros(2*n_pix + 1) + 1e-99

		### Loop over the evaluation points
		for jj in range(0,2*n_pix + 1):
			
			### Find the co-ordinates of the pixel in which the evaluation point lies
			xn = int(x_norm[jj])
			yn = int(y_norm[jj])

			### Loop over the pixels near the evaluation point
			for xx in range(xn-eval_dist,xn+eval_dist+1):
				for yy in range(yn-eval_dist,yn+eval_dist+1):

					### Calculate distance between the evaluation point and the pixels nearby, and therefore the pixel's weight
					dist = numpy.sqrt( (xx-x_norm[jj])**2 + (yy-y_norm[jj])**2 )
					weight = (1.0/numpy.sqrt(2*numpy.pi*h_length**2)) * numpy.exp(-(dist**2) / (2*h_length**2))

					### Add to the profile and the weight
					profile[jj] = profile[jj] + Map[xx,yy]*weight
					profile_weight[jj] = profile_weight[jj] + weight


		### Take the weighted mean of the contributions along the profile 
		profile = profile / profile_weight

		### Flip the profile if the gradient of the normal is negative
		if(-a/b > 0):
			profile = profile[::-1]
			
		Str_fil[:,ii] = profile


	### Construct radial position array
	radial = numpy.linspace(-max_dist,max_dist,2*n_pix+1)

	return Str_fil, length, radial










### A function to straighten a filament using interpolation
### Inputs:
### spine        (2d-array)             An ordered list of x,y co-ordinates of the spine points, shape n by 2 for n spine points
### Map          (2d-array)             The column density or integrated intensity map of the filament which you wish to straighten
### n_pix        (int)                  The number of evaluation points on either side of the radial profile
### max_dist     (float)                The distance to which the radial profile extends
### order        (int)                  The order of polynomial used to fit the spine points      
###
### Str_fil      (2d-array)             The straightened filament map
### length       (array)                The length of the straightened filament
### radial       (array)                The radius at each evaluation point

def Straighten_filament_interp(spine, Map, n_pix, max_dist, order): 

	### Unpack the spine array and determine the number of spine points
	x = spine[:,0]
	y = spine[:,1]
	n_spine = len(x)

	### Fit a polynomial to the spine points, then find the gradient by differentiating the function 
	p = numpy.polyfit(x,y,order)
	f = numpy.poly1d(p)
	p2 = p[:-1]
	q = numpy.arange(1,order+1,1)
	grad = numpy.poly1d(p2*q[::-1])

	### Create an empty array for the straightened filament
	Str_fil = numpy.zeros((2*n_pix + 1,n_spine))

	### Create empty length array
	length = numpy.zeros(n_spine)

	### Loop over the spine points
	for ii in range(0,n_spine):

		### If the spine point is not the first then we construct the length along the filament
		if(ii>0):
			length[ii] = length[ii-1] + numpy.sqrt( (x[ii]-x[ii-1])**2 + (y[ii]-y[ii-1])**2 )

		### Store co-ordinates of the spine point
		x0 = int(x[ii])
		y0 = int(y[ii])

		### Set up profile
		profile = numpy.zeros(2*n_pix + 1)
	
		### Find the gradient at this spine point and using it as a normal define a plane perpendicular to the spine
		a = 1
		b = grad(x0)
		norm = numpy.array([a,b],dtype=numpy.double)
		plane_const = numpy.sum(norm*spine[ii,:]) 

		### Construct the evaluation points along the normal line
		m = -a/b
		xmax = x0 + numpy.sqrt(max_dist**2/(m*m + 1))
		xmin = x0 - numpy.sqrt(max_dist**2/(m*m + 1))
		x_norm = numpy.linspace(xmin,xmax,2*n_pix + 1)
		y_norm = m*x_norm - m*x0 + y0

		### Loop over the evaluation points
		for jj in range(0,2*n_pix + 1):
			
			### Find the co-ordinates of the pixel in which the evaluation point lies
			xx = int(x_norm[jj])
			yy = int(y_norm[jj])

			### Calculate the derivatives for the interpolation
			dfdx = (Map[xx+1,yy] - Map[xx-1,yy])/2
			d2fdx2 = (Map[xx+2,yy] + Map[xx-2,yy] - 2*Map[xx,yy])/4

			dfdy = (Map[xx,yy+1] - Map[xx,yy-1])/2
			d2fdy2 = (Map[xx,yy+2] + Map[xx,yy-2] - 2*Map[xx,yy])/4

			d2fdxdy = (Map[xx+1,yy+1] - Map[xx-1,yy+1] - Map[xx+1,yy-1] + Map[xx-1,yy-1])/4

			### Find the distance between the pixel and the evaluation point
			dx = x_norm[jj] - xx
			dy = y_norm[jj] - yy

			### Second order Taylor expansion
			profile[jj] = Map[xx,yy] + dfdx*dx + dfdy*dy + 0.5*(d2fdx2*dx**2 + 2*d2fdxdy*dx*dy + d2fdy2*dy**2)

		### Flip the profile if the gradient of the normal is negative
		if(-a/b > 0):
			profile = profile[::-1]
			
		Str_fil[:,ii] = profile


	### Construct radial position array
	radial = numpy.linspace(-max_dist,max_dist,2*n_pix+1)

	return Str_fil, length, radial










### A function map core position onto the straightened filament
### Inputs:
### spine        (2d-array)             An ordered list of x,y co-ordinates of the spine points, shape n by 2 for n spine points
### pos          (2d-array)             The 2d array containing the x and y positions of cores associated with this filament     
### order        (int)                  The order of polynomial used to fit the spine points 
###
### core_pos     (2d-array)             The 2d array containing the x and y position of cores now mapped into the r-l co-ordinates of the straighten filament

def Map_cores(spine, pos, order): 

	### Make an empty array for the new core positions
	core_pos = numpy.zeros_like(pos)	
	n_cores = len(pos[:,0])

	### Extract spine positions and fit with a polynomial
	xs = spine[:,0]
	ys = spine[:,1]
	p = numpy.polyfit(xs,ys,order)
	f = numpy.poly1d(p)
	p2 = p[:-1]
	q = numpy.arange(1,order+1,1)
	grad = numpy.poly1d(p2*q[::-1])

	### Loop over each core and map
	for ii in range(0,n_cores):

		### Unpack x and y
		x = pos[ii,0]
		y = pos[ii,1]

		### Set min_dist to a large number and reset the index
		min_dist = 1e99
		min_dist_index = 1

		### Loop over all spine points to find the closest one
		for jj in range(0,len(spine[:,0])-1):
			
			### Unpack spine x and y
			xs = spine[jj,0]
			ys = spine[jj,1]

			### Determine distance between core and spine point
			dist = numpy.sqrt( (x-xs)**2 + (y-ys)**2 )

			### If distance is smaller than min distance then store it
			if(dist<min_dist):

				### Determine if the core appears on the right or the left and store this as the direction variable
				a = 1
				b = grad(xs)

				dxc = x-xs
				dyc = y-ys

				if(a*dyc - dxc*b > 0):
					direction = -1
				else:
					direction = 1

				min_dist = dist
				min_dist_index = jj

		### Once we go over all spine points we know the r-l co-ordinates of the core
		core_pos[ii,0]=direction * min_dist
		core_pos[ii,1]=min_dist_index
		
	return core_pos










########################################################################################
################# Function which are not called from outside FragMent ##################
########################################################################################



### A function to return the normalised PDF for two Gaussians
### Inputs:
### x             (array)               The evaluation points of the PDF
### params        (array)               The parameters describing the PDF
###
### Outputs:
### PDF           (array)               The PDF

def norm_two_guass(x, params):

	### Unpacking the parameters
	mu1, sig1, mu2, sig2, ratio = params

	### Determining the amplitudes of the two Gaussians 
	A1 = 1.0/(numpy.sqrt(2*numpy.pi*sig1*sig1) + ratio*numpy.sqrt(2*numpy.pi*sig2*sig2) )
	A2 = A1 * ratio

	### Calculating the PDF
	PDF = A1*numpy.exp(-(x - mu1)**2 / (2*sig1*sig1)) + A2*numpy.exp( -(x-mu2)**2 / (2*sig2*sig2))

	### Returning the PDF
	return PDF










### A function to return the log-likelihood of the data given a single Gaussian model
### Inputs:
### params        (array)               The parameters describing the model
### data          (array)               The data array
### bound         (array)               The boundaries for the model parameters
###
### Outputs:
### lnlike        (float)               The log-likelihood given that data and those model parameters

def lnlike_one_gauss(params, data, bound):

	### Unpack the parameters
	mean, sigma = params

	### Determine the number of data points
	n_data = len(data)

	### If the model parameters lie outside the boundaries for the parameters then return negative infinity
	if(mean < bound[0] or mean > bound[1] or sigma < bound[2] or sigma > bound[3]):
		return -numpy.inf

	### If model parameters lie within the boundaries then determine the log-likelihood
	else:
		lnlike = n_data * numpy.log(1.0/numpy.sqrt(2*numpy.pi*sigma*sigma)) - numpy.sum((data-mean)**2)/(2*sigma*sigma)
		return lnlike










### A function to return the log-likelihood of the data given the two Gaussian model
### Inputs:
### params        (array)               The parameters describing the model
### data          (array)               The data array
### bound         (array)               The boundaries for the model parameters
###
### Outputs:
### lnlike        (float)               The log-likelihood given that data and those model parameters

def lnlike_two_gauss(params, data, bound):

	### Unpack the parameters
	mu1, sig1, mu2, sig2, ratio = params

	### Determine the number of data points
	n_data = len(data)

	### If the model parameters lie outside the boundaries for the parameters then return negative infinity
	if(mu1 < bound[0] or mu1 > bound[1] or sig1 < bound[2] or sig1 > bound[3] or mu2 < bound[4] or mu2 > bound[5] or sig2 < bound[6] or sig2 > bound[7] or ratio < bound[8] or ratio > bound[9]):
		return -numpy.inf

	### If model parameters lie within the boundaries then determine the log-likelihood
	else:

		### Determine the amplitude of the two Gaussian's given model parameters
		A1 = 1.0/(numpy.sqrt(2*numpy.pi*sig1*sig1) + ratio*numpy.sqrt(2*numpy.pi*sig2*sig2) )
		A2 = A1 * ratio

		### Determine the likelihood of each data point and ensure none are zero
		like = A1*numpy.exp(-(data - mu1)**2 / (2*sig1*sig1)) + A2*numpy.exp( -(data-mu2)**2 / (2*sig2*sig2))
		like = like + 1e-99

		### Log the likelihood and sum to find the log-likelihood of the whole data-set
		lnlike = numpy.log(like)
		lnlike = numpy.sum(lnlike)

		return lnlike










### A function to return the un-normalised posterior of the data given the single Gaussian model
### Inputs:
### params        (array)               The parameters describing the model
### data          (array)               The data array
### bound         (array)               The boundaries for the model parameters
### prior         (function)            The prior for the model parameters
###
### Outputs:
### post          (float)               The un-normalised posterior for the given data-set and prior 

def post_one_gauss(params, data, bound, prior):

	### Unpack the parameters
	mu1, sig1 = params

	### If the model parameters lie outside the boundaries for the parameters then return negative infinity
	if(mu1 < bound[0] or mu1 > bound[1] or sig1 < bound[2] or sig1 > bound[3]):
		return -numpy.inf

	### If model parameters lie within the boundaries then determine the posterior
	else:
		### Determine value of prior and data-likelihood for these parameters
		p = prior(params, bound)
		l = lnlike_one_gauss(params, data, bound)

		### Combine to work out un-normalised posterior
		post = p * numpy.exp(l)

		### Return value
		return post










### A function to return the un-normalised posterior of the data given the two Gaussian model
### Inputs:
### params        (array)               The parameters describing the model
### data          (array)               The data array
### bound         (array)               The boundaries for the model parameters
### prior         (function)            The prior for the model parameters
###
### Outputs:
### post          (float)               The un-normalised posterior for the given data-set and prior 

def post_two_gauss(params, data, bound, prior):

	### Unpack the parameters
	mu1, sig1, mu2, sig2, ratio = params

	### If the model parameters lie outside the boundaries for the parameters then return negative infinity
	if(mu1 < bound[0] or mu1 > bound[1] or sig1 < bound[2] or sig1 > bound[3] or mu2 < bound[4] or mu2 > bound[5] or sig2 < bound[6] or sig2 > bound[7] or ratio < bound[8] or ratio > bound[9]):
		return -numpy.inf

	### If model parameters lie within the boundaries then determine the posterior
	else:
		### Determine value of prior and data-likelihood for these parameters
		p = prior(params, bound)
		l = lnlike_two_gauss(params, data, bound)

		### Combine to work out un-normalised posterior
		post = p * numpy.exp(l)

		### Return value
		return post
