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

```
gcc -shared -Wl,-soname,FragMent_C -o FragMent_C.so -fPIC FragMent_C.c -g
```
