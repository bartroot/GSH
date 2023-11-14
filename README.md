# Global Spherical Harmonic package (GSH) 

GSH is a MATLAB package to do **Global Spherical Harmonic Analyses** (GSHA) and **Synthesis** (GSHS) for Crust1.0.


## Requirements

The code runs from Matlab 2016a, but highest version is recommended. 

SHBundle (https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/) is not needed to be downloaded, however the GSH code uses some functions from this bundle:

- cs2sc.m
- sc2cs.m
- gsha.m (has been altered for speed and clarity of the code by Bart Root)

This SHbundle is licenced under NU General Public Licence: (https://www.gnu.org/licenses/gpl-3.0.en.html). We highly recommend to download this bundle as well, because it has many usefull Matlab functions that can be used together with the GSH package.

## Structure

In this repository you can find:

- `run.m` executable to run the **GSHA** module (performed by `model_SH_analysis.m` found in `Tools/`) and **GSHS** module (performed by `model_SH_synthesis.m` found in `Tools/`)  
- `inputModel.m` used as input by the GSHA module   
- `plotResults_full.m` executable that can be used to visualize the results  
- `Data/` directory containing GMT files with the density model specifications   
- `Tools/` directory containing the scripts used to perform the GSHA and GSHS  
- `Results/` directory containing the outputs of `run.m`  

`run.m` can run the GSHA and the GSHS modules separately. THE GSHA module uses a `importfile = inputModel.m`. Here the density model can be defined by using .gmt files. GMT files must have a certain format, which can be inspected in the given tutorial files (see .gmt files in `Data/`). The unit in .gmt files are SI/1000 (km or g/cm<sup>3</sup>).

The GSHS module will return a structure `data` containing the results. These results can be visualized with `plotResults_full.m`

The `run.m` file can be modified if different area or resolution is wanted.

`run.m` uses: 

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


### Extra Tools

Other scripts that are found in `Tools/`:

- `visual_gmtfile.m`: reads a .gmt file and plots the file. It allows you to check the file and, if an error is present in the file, it will not run (then the GSH will also not work).   
- `degreeVariance.m`: generates the degree variances from the V vector, containing the spherical harmonic coefficients from the analysis.  
- `Europe_centered.m`: transforms a map matrix and swaps from -180:180 to 0:360 longitude view, as different conventions are available.  
- `matrix2gmt.m` and `gmt2matrix.m`: convert matrixes to .gmt structure and vice versa respectively.
- `topo2crust.m`: creates a flexure type crust-mantle interface from topography map (Airy, infinite plate, thin shell)
- `segment_2_layer_model.m`: makes segments of a 2 layer model to counter numerical instabilities.
- `GSHS.m`: Counterfunction of GSHA, synthesis coefficients to map on a unit sphere
- `sc2zdenek.m`: conversion of Stokes coefficients from GSH format to SFEC format
- `zdenek2V.m`: conversion of Stokes coefficients from SFEC format to GSH format


### Keep in Mind

Common error is a 'single' value `NaN` in the .gmt file, that arises due to interpolation error. This needs to be checked, because the GSH code will also stop working, when this is the case. Check all input files and variables for NaN values. This could happen if interpolation of the input maps has gone wrong. The isnan.m function of Matlab is very useful for this.

## Authors (Standing on the shoulders of giants!)

This Software has been developed on ideas and software from the following developers/contributors:

- **Bart Root** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0001-7742-1434](https://orcid.org/0000-0001-7742-1434), Technische Universiteit Delft (Developer)
- **Pavel Novak**, University of West Bohemia (Contributor)  
- **Dimitri Tsoulis**, Institut für Astronomische und Physikalische Geodäsie (IAPG), Technische Universität München  (Contributor)  
- **Nico Sneeuw** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0003-1796-0131](https://orcid.org/0000-0003-1796-0131), Institut für Astronomische und Physikalische Geodäsie (IAPG), Technische Universität München and University of Stuttgart (Contributor)  
- **Matthias Weigelt** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0001-9669-127X](https://orcid.org/0000-0001-9669-127X), Institut für Astronomische und Physikalische Geodäsie (IAPG), Technische Universität München  (Contributor)
- **Wouter van der Wal** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0001-8030-9080](https://orcid.org/0000-0001-8030-9080), Technische Universiteit Delft (Contributor)
- **Weilun Qin** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0001-7703-3132](https://orcid.org/0000-0001-7703-3132), Technische Universiteit Delft (Contributor)


## License

The contents of this repository are licensed under a **GNU General Public License v3.0** (see [LICENSE](https://github.com/bartroot/GSH/blob/main/LICENSE.md) file).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program “GSH package”. GSH is a MATLAB package to do Global Spherical Harmonic Analyses (GSHA) and Synthesis (GSHS) for Crust1.0. written by the Author(s). 
Henri Werij, Faculty of Aerospace Engineering, Technische Universiteit Delft. 

© 2021, B.C. Root

## References

If you would like to reuse the code, please cite the following DOI: 10.4121/16764238

The mathematical description can be found in the following publication:

Root, B.C., Novák, P., Dirkx, D., Kaban, M., van der Wal, W., and Vermeersen, L.L.A. (2016), On a spectral method for forward gravity field modelling, Journal of Geodynamics, 97, 22-30 [DOI: 10.1016/j.jog.2016.02.008](https://doi.org/10.1016/j.jog.2016.02.008)


## Would you like to contribute?

If you have any comments, feedback, or recommendations, feel free to **reach out** by sending an email to b.c.root@tudelft.nl

If you would like to contribute directly, you are welcome to **fork** this repository.

Thank you and enjoy!

