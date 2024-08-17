# photGalIMF version 1.2

Authors: Moritz Haslbauer (Bonn University), Yan Zhiqiang (Nanjing University)

The photGalIMF code can calculate the stellar mass and luminosity evolution of a galaxy with known information on age, mass, metallicity, and IMF of all the single stellar populations ever formed. 

That is, this is a post-processing code. It is a supplementary repository of the GalIMF galaxy chemical evolution code (https://github.com/Azeret/galIMF) to calculate (post-processing) the photometric luminosity of simulated galaxies.

## Stellar isochrone

Currently, the only stellar isochrone prepared is from the PARSEC stellar evolution model downloaded from http://stev.oapd.inaf.it/cgi-bin/cmd with settings detailed in the publication below.

## Citation

Haslbauer et al. (2024 A&A accepted)

## Input

The photGalIMF code takes the galaxy chemical evolution simulation results saved in (https://github.com/Azeret/galIMF/tree/master/simulation_results_from_galaxy_evol) and the galaxy-wide IMFs saved in (https://github.com/Azeret/galIMF/tree/master/Generated_IGIMFs) as inputs.

## Output

The photGalIMF code outputs the stellar mass evolution, dynamical mass (stellar + remnant mass) evolution, total stellar luminosity evolution (for V, Ks or IRAC36 bands), and stellar mass-to-light ratio evolution (for V, Ks or IRAC36 bands).

## Caveat

For the current luminosity table (to be updated), the minimum time resolution is 10 Myr and the minimum stellar population age is also 10 Myr. Therefore, the results are inaccurate when there are recent star formation activities.

## Example

This repository includes an example of galaxy simulation results and is self-sufficient to run without the GalIMF code.

You could test the code by directly running
stellar_luminosity.py
and/or
galaxy_luminosity_evolution.py

When a new galaxy simulation rather than the example is considered. Follow these steps to modify "galaxy_luminosity_evolution.py".

### Step 1

Specify the galaxy-wide IMF model (parameter "gwIMF") for the function "mass_to_light". 
if gwIMF == "Kroupa": The code adopts the invariant Kroupa 2001 IMF given in Kroupa_IMF.py
Else: The code adopts variable gwIMF for each timestep defined in simulation_results_from_galaxy_evol/igimf_epoch_{}.py

### Step 2

For variable gwIMF (see Step 1) the gwIMF for each timestep is required in simulation_results_from_galaxy_evol (see Inputs above).

### Step 3

Modify the file path. That is, the following lines:
```
sys.path.append('simulation_results_from_galaxy_evol/example')
file_evo = open('simulation_results_from_galaxy_evol/example/chemical_and_SN_evolution.txt', 'r')
```

### Step 4

Specify the following lines in "stellar_luminosity.py"
```
isochrones_selection = "PARSEC"
band_selection = "V"  # can be "V", "Ks", or "IRAC36"
```

### Step 5

Run galaxy_luminosity_evolution.py

## Updates

14.8.2024: 
The previous line index for reading Z__list, simulation_time_list, stellar_mass_formed_at_each_epoch are incorrect. Line index finders are implemented in galaxy_luminosity_evolution.py to automatically prepare these lists and also epoch_number_list so that it does not need to be adjusted manually anymore.

16.8.2024: 
Wrong luminosity was returned by function SSP_luminosity in stellar_luminosity.py, when the stellar metallicity is exactly one of the values in the Z_list_value list. This is now corrected.
This update does not affect calculations for galaxies with natural self-enrichment (and Haslbauer et al. 2024) because the metallicities given by a GCE model can rarely be exactly a number in the Z_list_value list.

