# Plotting Absorption Spectrum

This code plots the absorption spectrum from the output file (EXC.DAT) generated by the DFTB+ 1.3 code.

## Running DFTB+ calculations

A sample input for the DFTB+ code and the coordinate file is given for the benzene molecule. After the calculation completes, the DFTB+ code will create the EXC.DAT file

This file contains in columns, the excitation energies and the corresponding oscillator strengths.

## Running the python code

Run the python code once you have the EXC.DAT file. The code will show a plot of the absorption spectrum and also save the information to a text file, if you prefer some other code to plot the spectrum
