# photGalIMF

## Aim

Supplementary repository of the GalIMF galaxy chemical evolution code (https://github.com/Azeret/galIMF) to calculate (post-processing) the photometric luminosity of simulated galaxies.

## Documentation

This is kept as an independent repository that may be updated with newer stellar evolution models.

## Input

The photGalIMF code takes the galaxy chemical evolution simulation results saved in (https://github.com/Azeret/galIMF/tree/master/simulation_results_from_galaxy_evol) and the galaxy-wide IMFs saved in (https://github.com/Azeret/galIMF/tree/master/Generated_IGIMFs) as inputs.

## Output

The photGalIMF code outputs the stellar mass evolution, dynamical mass (stellar + remnant mass) evolution, total stellar luminosity evolution (for V, Ks or IRAC36 bands), and stellar mass-to-light ratio evolution (for V, Ks or IRAC36 bands).

## Caveat

For the current luminosity table (to be updated), the minimum time resolution is 10 Myr and the minimum stellar population age is also 10 Myr. Therefore, the results are inaccurate when there are recent star formation activities.

## Example

This repository includes an example of galaxy simulation results and is self-sufficient to run without GalIMF.

You could test the code by directly running
stellar_luminosity.py
and/or
galaxy_luminosity_evolution.py

If a new galaxy simulation rather than the example is considered. Follow these steps to modify "galaxy_luminosity_evolution.py".

### Step 1

Specify the galaxy-wide IMF model (parameter "gwIMF") for the function "mass_to_light". 
if gwIMF == "Kroupa": The code adopts the invariant Kroupa 2001 IMF given in Kroupa_IMF.py
Else: The code adopts variable gwIMF for each timestep defined in simulation_results_from_galaxy_evol/igimf_epoch_{}.py

### Step 2

Specify the "epoch_number_list", which is a list of timesteps (real-time / 10 Myr) that have star formation activities. 
The corresponding total stellar mass formed and their mean stellar metallicity for each timestep is provided in simulation_results_from_galaxy_evol/chemical_and_SN_evolution.txt (see Inputs above).
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
