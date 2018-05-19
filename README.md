Documentation for Github data repository:

================================================================
 
 Code name: main.py
 
 Code author: Khai Nguyen
 
 Date: 18 May 2018
 
 Language: Python
 
================================================================

This Python script is a resource accompanying the following manuscript: Bryan J. Pflueger, Khai Nguyen, Tamara Bogdanovic, Michael Eracleous, Jessie C. Runnoe, Steinn Sigurdsson, and Todd Boroson 2018, “Likelihood for Detection of sub-parsec SBHBs in Spectroscopic Surveys”, ApJ, submitted (arXiv:1803.02368). It is provided as is, free of charge and with no technical support. We ask you to please cite the original manuscript if you use the script or any of the manuscript data.

The script produces a 2D likelihood map as a function of the mass ratio and orbital semi-major axis, as shown in Figure 8 of Pflueger et al. (2018). The figure is created as a PDF file with a default name <Likelihood_q_a.pdf>. The script includes three modules, the main module <main.py>, the module called by the main module <functions.py>, as well as the input ascii file <input.dat>, which contains the model parameters and constants used in the calculation.

The header of the <input.dat> file includes the description of units, format, and nominal values for the relevant parameters. Note that binary mass is defined as an array of values (mlist) and so is the mass ratio (mdotlist). The likelihood figure, which is calculated for a given binary mass and mass accretion rate, is plotted for the first value in the mlist and mdotlist arrays.

