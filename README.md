picarro
=======

R script for processing output files from a Picarro G2301.

You can point this script to a directory of output files downloaded from the G2301, an infrared CH4/CO2/H2O analyzer made by Picarro. It will process them: identifying outliers, computing gas fluxes, summarizing files, and writing outputs and logs.

The script is written for R 3.0.3 and doesn't use anything particularly exotic. Required packages are listed in the script.

Ben Bond-Lamberty
bondlamberty@pnnl.gov
