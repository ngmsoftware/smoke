# Smoke

Matlab smoke simulator based on numerical solution of Navier-Stokes equations. This is just for entertainment purposes. The numerical solution is not precise and involves "empirical" steps like smoothing of the fields, etc.
Also, its sensible to instabilization if the parameters are not well behaved

## Software

Matlab

## Hardware

## Screen-shoots

![Screenshoot 1](/doc/img1.png?raw=true "Sample image 1")
![Screenshoot 2](/doc/img2.png?raw=true "Sample image 2")
![Screenshoot 3](/doc/img3.png?raw=true "Sample image 3")

## Instructions

1 - If necessary, compile the navierStokesStep.c file using 

\>\> mex navierStokesStep.c

2 - Adjust the parameters (source position, intensity, etc...) in the SOURCE and FIELD configuration section 

3 - execute run.m