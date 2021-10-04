# GSH
Global Spherical Harmonic package (MATLAB)

README: Gravity Software TUDelft

NOTE: This software is not tested thoroughly for commercial use. If used in research please give credits to the developer and his University. For any other License questions, please consult the LICENSE.md file.

Initial development by Bart Root, 06-nov-2014, Delft University of Technology:

The following software package contains MATLAB scripts to do Global Spherical  Harmonic Analyses and Synthesis for Crust1.0. The setup is build as follows:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
There are two executables

run.m

and 

plotResults_full.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run.m is divided in to parts:

1. GSHA -> performed by model_SH_analysis.m
2. GSHS -> performed by model_SH_synthesis.m

Both modules can be used separate. THE GSHA module uses a importfile = inputModel.m. Here the density model can be defined by using .gmt files. GMT files must have a certain format, which can be inspected in the given tutorial files (dir: Data). The unit in .gmt files are SI/1000. So, km or g/cm3.

The GSHS module will return a structure ‘data’. In this structure are all the different results situated. These results can be view with plotResults_full.m

The run.m file can be modified if different area or resolution is wanted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Structure:

Run.m uses: 

	- inputModel.m
	- model_SH_synthesis.m
		- gravityModule.m
			- Legendre_functions.m
		- gravityModule_full.m (basic tested, slow!!!!)
			- Legendre_functions.m
	- model_SH_analysis.m
		- import_layer.m
			- gmt2matrix.m
		- layer_SH_analysis.m
			- GSHA.m
				-cs2sc.m
				-Legendre_functions.m
			- sc2vecml.m
			- geocradius.m (MATLAB func. toolbox aero)
					and only when geoid is WGS84

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EXTRA TOOLS

Another useful script is the visual_gmtfile.m located in the Tools directory.
This function will read in a gmt file and plot the file. It allows you to check 
easily the file and if an error is present in the file it will not run, but then the 
GSH code will also not work. 

degreeVariance.m generates the degree variances from the V vector, containing the spherical harmonic coefficients from the analysis.

Europe_centered.m transforms a map matrix and swaps from -180:180 to 0:360 longitude view, as different conventions are available.

gmt2matrix.m and matrix2gmt.m convert matrixes to .gmt structure and visa versa.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Common error is a 'single' value NaN in the .gmt file, due to interpolation error.
This needs to be checked, because the GSH code will also stop working, when this is the case.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTRIBUTIONS (Standing on the shoulders of giants)

This Software has been developed on ideas and software from the following developers/contributors:

- Pavel Novak, University of West Bohemia (Developer)
- Dimitri Tsoulis, IAPG, TU-Munich (Developer)
- Nico Sneeuw, IAPG, TU-Munich (Developer): https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/ 
- Matthias Weigelt, IAPG, TU-Munich (Developer)
- Wouter van der Wal, Delft University of Technology (Contributor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REFERENCES

This code is based on the mathematical representation of:

Root, B.C., Novák, P., Dirkx, D., Kaban, M., van der Wal, W., and Vermeersen, L.L.A. (2016), On a spectral method for forward gravity field modelling, Journal of Geodynamics, 97, 22-30.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

If you have any comments, feedback, or recommendations, please let me know: b.c.root@tudelft.nl

Thank you and enjoy!
