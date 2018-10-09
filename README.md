# A code to study the fragmentation of filaments (FragMent)

FragMent is a python/C module which collates a number techniques used to study fragmentation in filaments. It also performs model selection using a frequentist and Bayesian approach to find the best descriptor of a filament's fragmentation. The code is open-source and can be downloaded here. 

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

### Nearest neighbour separations

The nearest neighbour method is called by using the function
```
nn_sep = Fragment.NNS(pos)
```
where *pos* is the array containing the x and y core positions, and *nn_sep* is the returned 1D array containing the nearest neighbour separations. The results are then taken as inputs for 3 null hypothesis tests:
```
p_median, p_mean = FragMent.Stat_Sig(nn_seps, "NNS", boundary, nruns, sep_limit)
ks_stat, p_ks = FragMent.KS_test(nn_seps, "NNS", boundary, nruns, sep_limit)
ad_stat, crit_vals, p_ad = FragMent.AD_test(nn_seps, "NNS", boundary, nruns, sep_limit)
```
The option "NNS" is included to tell each function that the nearest neighbour separation method is being tested. The *boundary* array is as described above. *nruns* is the number of random realisations used to construct the null hypothesis, a typical value is 10,000, although the higher the better. *sep_limit* is the smallest separation that can be detected; for observations it would be 1 or 2 beamsizes. These functions return the p-values for each of the tests. For the Kolmogorov-Smirnov test, it also returns the K-S statistic. For the Anderson-Darling test it returns the A-D statistic as well as a critical value array which shows the values of the A-D statistic corresponding to the significance levels: 25%, 10%, 5%, 2.5%, 1%.

### Minimum spanning tree

The minimum spanning tree method is used similarly to the nearest neighbour method,
```
mst_seps, mst = FragMent.MST(pos)
```
where *mst_seps* is the 1D array containing the edge lengths, and *mst* is the 2D array containing the minimum spanning tree.

The same three functions are used to perform null hypothesis tests but "MST" is used instead of "NNS".

### Two-point correlation function

The two-point correlation function is called using either an exact KDE
```
sep, w = FragMent.TwoPoint(pos,nruns,boundary,sep_limit)
```
or an approximate KDE
'''
sep, w = FragMent.ApproxTwoPoint(pos,nruns,boundary,sep_limit,op_crit)
'''
The *pos*, *nruns*, *boundary* and *sep_limit* are as before. *op_crit* is the critical opening distance used for the tree when constructing the approximate KDEs. Its value depends on the data size, *op_crit* > 10<sup>4</sup> / *N*, where *N* is the number of data points used to construct the KDE, *N* = *nruns* x *ncores* <sup>2</sup>.



