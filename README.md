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
gcc -shared -Wl,-install_name,FragMent_C -o FragMent_C.so -fPIC FragMent_C.c -g
```

## Using the code

Included in the repository is an example script which calls all the functions in FragMent. Here I go through each function and detail their parameters.

### Required information 

Three pieces of information are needed for the five fragmentation techniques: the core positions, the map values along the spine, and the boundary of the map. 

* The core positions should be in a plain text file consisting on two columns, the x and y positions. All separations returned by the functions will be given in the same units as the core positions (i.e. degrees, arcseconds, parsec). These will be used by the nearest neighbour separation, minimum spanning tree, N<sup>th</sup> nearest neighbour, and two-point correlation function methods.

* The map values (column density / integrated intensity) along the spine should be in a plain text file. This is used for the Fourier transform function.

* The map boundary defines the boundary within which randomly placed cores are placed. It is defined by the maximum and minimum value in the x and y direction. 

### Nearest neighbour separations


