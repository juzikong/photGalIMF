# photGalIMF

Supplementary repository of the GalIMF galaxy chemical evolution code (https://github.com/Azeret/galIMF) to calculate the photometric luminosity of simulated galaxies.

This is kept as an independent repository that may be updated with newer stellar evolution models.

The photGalIMF code takes the galaxy chemical evolution simulation results saved in (https://github.com/Azeret/galIMF/tree/master/simulation_results_from_galaxy_evol) and the galaxy-wide IMFs saved in (https://github.com/Azeret/galIMF/tree/master/Generated_IGIMFs) as inputs.

The photGalIMF code outputs the luminosity and mass-to-light ratio evolutions.

This repository includes an example of galaxy simulation result and is self-sufficient to run without GalIMF.

You could test the code by directly running
stellar_luminosity.py
and/or
galaxy_luminosity_evolution.py
