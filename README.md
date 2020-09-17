![alt text](https://github.com/SeamusClarke/FragMent/blob/master/Images/FragMent.png)

# A code to study the fragmentation of filaments (FragMent)

FragMent is a python/C module which collates a number techniques used to study fragmentation in filaments. It also performs model selection using a frequentist and Bayesian approach to find the best descriptor of a filament's fragmentation. The code is open-source and can be downloaded here.

The accompanying paper can be found [here](http://adsabs.harvard.edu/abs/2019arXiv190106205C). It details the sensitivities of each method and explains in more detail a number of the procedures one should use when analysing fragmentation.

While the code was designed to investigate filament fragmentation the functions are general and may be used for any set of 2D points to study more general cases of fragmentation. 

## Dependencies 

FragMent requires 4 common libraries to be installed:

* Numpy,
* Scipy,
* ctypes,
* emcee. 

To allow the import of the FragMent module from any directory use the export command to modified the PYTHONPATH variable. This is done by adding the line

```
export PYTHONPATH=$PYTHONPATH:"Path to FragMent's download location"
```
to the .bashrc file in the home directory.

To compile the C portion of the module use the command

For Linux:
```
gcc -shared -Wl,-soname,FragMent_C -o FragMent_C.so -fPIC FragMent_C.c -g
```
For Mac:
```
cc -shared -Wl,-install_name,FragMent_C -o FragMent_C.so -fPIC FragMent_C.c -g
```

## Using the code

Included in the repository is an example script which calls all the functions in FragMent. Here I go through each method included in FragMent and detail how to use them.

### Required information 

Three pieces of information are needed for the five fragmentation techniques: the core positions, the map values along the spine, and the boundary of the map. 

* The core positions should be in a plain text file consisting on two columns, the x and y positions. All separations returned by the functions will be given in the same units as the core positions (i.e. degrees, arcseconds, parsec). These will be used by the nearest neighbour separation, minimum spanning tree, N<sup>th</sup> nearest neighbour, and two-point correlation function methods.

* The map values (column density / integrated intensity) along the spine should be in a plain text file. This is used for the Fourier transform function.

* The map boundary defines the boundary within which randomly placed cores are placed. It is set by defining maximum and minimum values in the x and y directions. 

### Fragmentation methods

#### Nearest neighbour separations

The nearest neighbour method is called by using the function
```
nn_sep = Fragment.NNS(pos)
```
where *pos* is the array containing the x and y core positions, and *nn_sep* is the returned 1D array containing the nearest neighbour separations.

#### Minimum spanning tree

The minimum spanning tree method is used similarly to the nearest neighbour method,
```
mst_seps, mst = FragMent.MST(pos)
```
where *mst_seps* is the 1D array containing the edge lengths, and *mst* is the 2D array containing the minimum spanning tree.

#### Two-point correlation function

The two-point correlation function is called using either an exact KDE
```
sep, w = FragMent.TwoPoint(pos,nruns,boundary,sep_limit)
```
or an approximate KDE
```
sep, w = FragMent.ApproxTwoPoint(pos,nruns,boundary,sep_limit,op_crit)
```
*pos*, *nruns*, *boundary* and *sep_limit* are as before. *op_crit* is the critical opening distance used for the tree when constructing the approximate KDEs. Its value depends on the data size, *op_crit* > 10<sup>4</sup> / *N*, where *N* is the number of data points used to construct the KDE, *N* = *nruns* x *ncores* <sup>2</sup>. The approximate function should typically be used as it is 50 to 100 times faster than the exact. 

Both functions return the separation array, *sep*, the separations at which the two-point correlation function is evaluated; and *w*, the two-point correlation function itself.

#### N<sup>th</sup> nearest neighbour

The N<sup>th</sup> nearest neighbour is called using the function
```
NNN_sep = FragMent.NNNS(pos)
```
*NNN_sep* is a 2D array containing the neighbour separations for each core. It is reduced to the mean and standard deviation for each N<sup>th</sup> neighbour using the following lines
```
mean_seps = numpy.mean(NNN_sep, axis=0)
std_seps = numpy.std(NNN_sep, axis=0)
```

#### Fourier power spectrum

The Fourier power spectrum requires the map values along the spine. It is called by the command
```
k,power = FragMent.FT_Spine(spine)
```
It returns the array of wavenumbers, *k*, and the power spectrum at those wavenumbers, *power*. The first element of both of these arrays must be removed as it corresponds to the *k*=0 mode. The length scale, *l*, that corresponds to wavenumber *k* is here given by *l*=*L*/*k*, where *L* is the length of the spine. 

### Statistical significance tests

Three functions are used to perform null hypothesis tests of the results from the nearest neighbour and minimum spanning tree method. These construct distributions assuming the null that there is no correlation between core locations, i.e. they are random, and return p-values.

#### Kolmogorov-Smirnov test 

The Kolmogorov-Smirnov test is called using the function
```
stat_dist,p_val = FragMent.KS_test(dist,method,boundary,num,lower_lim)
```
It returns two variables: *stat_dist*, the K-S test statistic; and *p_val*, the p-value associated with the test statistic given the number of data points. It takes 5 inputs: *dist*, the array of separations generating either from the nearest neighbour of minimum spanning tree methods; *method*, a string set to "NNS" if the array of separations is from the nearest neighbour method or "MST" if from the minimum spanning tree method; *boundary*, the boundary array denoting the area in which to randomly place cores; *num*, the number of random realisations used to construct the null hypothesis, a typical value is 10,000, although the higher the better (but slower); *lower_lim*, is the smallest separation that can be detected, i.e. 1 or 2 beamsizes. 

#### Anderson-Darling test

The Anderson-Darling test is called similarly to the KS test: 
```
ad_stat, crit_vals, p_ad = FragMent.AD_test(dist,method,boundary,num,lower_lim)
```
It has the same 5 inputs as the KS test. It returns three variables: *ad_stat*, the Anderson-Darling test statistic; *crit_vals*, is an array containing the values of the A-D test statistic values corresponding to the significant levels 25%, 10%, 5%, 2.5% and 1%; and *p_ad*, the p-value of the null-hypothesis test.

#### Average and with null hypothesis test

As well as the KS and AD tests, FragMent includes a null hypothesis test based on the average and width measurements of the separation distributions, as described in the FragMent paper. This test is called using the command:
```
p_median, p_mean = FragMent.Stat_Sig(dist,method,boundary,num,lower_lim)
```
The function uses the same 5 inputs as the KS and AD test. It has two output: *p_median*, the p-value when the median and the interquartile range are used as the average and width estimators; and *p_mean*, the p-value when the mean and standard deviation are used. 

### Model selection techniques

If the null can be rejected the next step is to determine if the fragmentation can be best characterised by one or two fragmentation length-scales. FragMent includes a frequentist and a Bayesian approach to this problem.

#### Akaike information criterion - a frequentist approach

The Akaike information criterion (AIC) can used to select a best fitting model. This is done using the function:
```
model_one, theta_one, model_two, theta_two = FragMent.Freq_model_selection(mst_seps,bound,mcmc_param)
```
This function requires 3 inputs: *mst_seps*, the separation distribution from the minimum spanning tree method; *bound*, an array containing the boundaries for the five parameters of the two-tier model; and *mcmc_param*, which is an array containing the Monte Carlo Markov Chain parameters to be used.

The *bound* array is 10 entries long, describing the minimum and maximum values allowed for the 5 parameters of the two-tier fragmentation model: the mean of the first Gaussian, the standard deviation of the first Gaussian, the mean of the second Gaussian, the standard deviation of the second Gaussian, and the ratio of the amplitudes of the two Gaussians. 

The *mcmc_param* array is 3 entries long and of type int. The first entry sets the number of walkers to be used; the second entry the number of burn steps that each walker takes; and the third entry sets the number of steps the walkers take after the burn stage.

The output of the function are the two *model* arrays and the two *theta* arrays. The *model* arrays contain the AIC of the model as the first entry, and the Akaike weight as the second entry. The *theta* arrays contain the fitting model parameters which maximise the data likelihood. 

#### Odds ratio - a Bayesian approach

The Bayesian evidence of each model is currently calculated explicitly. This is done using the function:
```
evi_one, evi_two = FragMent.Evidence_explicit(mst_seps,bound,prior_one, prior_two,n)
```
The function takes 5 inputs. *mst_seps* and *bounds* are the same as those for the Akaike information criterion. *prior_one* and *prior_two* are functions which calculate the priors for model one and two given a set of parameters. *n* is the number of grid points in 1D which are used to calculate the posterior, *n* should be of the order 30 or 40. 

The *prior* functions must take two arguments: *params*, an array containing the current value of the model parameters; and *bounds*, the same array as that passed to *Evidence_explicit*.

The function returns 2 variables, the Bayesian evidence of each model. The odds ratio is then calculated by dividing one by the other.

### Straightening filaments

One should straighten curved filaments before applying the fragmentation techniques included in FragMent, this is to remove the complexity added by the curvature. FragMent includes two different ways to do this.

#### Interpolation

The first is to use interpolation to determine the radial profile at every spine point. This is done using the function:
```
Str_fil, length, radial = FragMent.Straighten_filament_interp(spine, Map, n_pix, max_dist, order)
```
The function returns 3 outputs: *Str_fil*, the map of the straightened filament; *length*, a one-dimensional array detailing the length of the filament at each spine pixel; *radial*, a one-dimensional array detailing the radial distance from the spine at each radial position.

The 5 inputs are: *spine*, the pixel co-ordinates of the filament spine given in order from one end of the spine to the other; *Map*, the column density or integrated intensity map; *n_pix*, the number of evaluation points on either side of the spine for the radial profile; *max_dist*, the maximum distance to which the radial profile extends, given in pixels; and *order*, the order of the polynomial used to fit the spine points, typically ~10. 

#### Kernel weighting

The second method uses a Gaussian kernel to determine the radial profile at each spine point. This is slower than the interpolation function but is more stable when there are large gradients in the map. It is called using the function:
```
Str_fil, length, radial = FragMent.Straighten_filament_weight(spine, Map, n_pix, max_dist, order, h_length)
```
It has the same inputs and outputs as the interpolation method but includes one more input, *h_length*. This is the smoothing length for the kernel, measured in pixels. A typical value is 0.5 pixels.

#### Core mapping

Once the filaments are straightened the cores must be mapped into this straightened filament space. This is done by using the function:
```
core_pos = FragMent.Map_cores(spine, pos, order) 
```
The inputs are: *spine*, the ordered list of spine pixel values used in the straightening function; *pos* the 2D array of core positions, and *order*, the order of the polynomial used to fit the spine points when straightening. The output is *core_pos*, the 2D array of core positions in the length-radius space of the straight filament. It is these positions which should act as inputs to the fragmentation analysis functions.

## Acknowledgements 
This code was produced with the support of the ERC starting grant No. 679852 "RADFEEDBACK"

## Contact

Dr. Seamus Clarke <br /> 
I. Physikalisches Institut, <br /> 
Universität zu Köln, <br /> 
Zülpicher Str. 77, <br /> 
D-50937 Köln, <br /> 
Germany

clarke at ph1.uni-koeln.de 


