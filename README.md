# cosmo-inversion
MCMC inversion code for cosmogenic nuclide applications
This MCMC inversion code generates the model output presented and discussed in the manuscript: A topographic hinge-zone divides coastal and inland ice dynamic regimes in East Antarctica, Andersen et al. 2023: https://doi.org/10.1038/s43247-022-00673-6 

Note 1: At present the code is set up to handle 10Be and 26Al (mandatory), while 36Cl and 21Ne are optional

Note 2: Production rates are presently calculated using the 'St' scaling scheme. This is primarily due to the lack of an open source framework for calculating the 'LSDn' production rates for all four nuclides in a computationally consistent way. 'LSDn' scaling is 
better, especially at high latitudes.

Primary functions
1. compile_mn_data.m: Reads metadata and multi-cosmogenic nuclide data for a set of samples and saves data structure including production parameters in .mat format.

2. bedrockMCvJ1_mn5.m: Performs inversion for one or more samples. If multiple samples are provided they are assumed to have the same exposure history while the exhumation history is allowed to vary between samples. To test performance reduce number of walkers/number of models per walker in the code. To invert samples used in main manuscript call this function in a for-loop (consider parallelising):

	for i=1:34, bedrockMCvJ1_mn5(i,2,0,1); end %Heimefrontfjella samples
	
	
	for i=1:29,bedrockMCvJ1_mn5(i,2,0,2); end %Jutulstraumen samples
	
	
 Output is saved to the folder models/MDML
	
3-6. forward_bedrockvJ1_mn4_BeAl_xx.m: Set of four m-files to forward calculate 		nuclide concentrations in a sample given the exposure and exhumation history 		parameters. If ismember lines are uncommented only one file is needed to calculate 	all four nuclides, but the code will be very slow.

7. compile_results_mn5.m: Compiles inversion model results for all samples from each site for use in plotting scripts.

8-9. MDMLx_elevation_fig.m: These two m-files read and plot compiled results for 		Heimefrontfjella (MDML1) and Jutulstraumen (MDML2) ordered by elevation above 	sea level.

10. makereportJ_mn5.m: Plots information about walkers, parameter distributions, model-data fit, and saves report as pdf in models/MDML/reports folder.

Contents of subfolders
- The data folder contains input_mn.xlsx and stores mat-files generated from compile_mn_data.m
- The ClimateCurves folder contains:
    - d18Ocurves.mat (Lisiecki & Raymo 2005), Data downloaded from https://lorraine-lisiecki.com/stack.html
    - zachos_.mat (Zachos et al., 2001; default in manuscript),
    - milankovitch.mat. Data described in Laskar 2011 (La2004) and downloaded here: https://biocycle.atmos.colostate.edu/shiny/Milankovitch/
    - iceDMC.mat, Ice-sheet model output data generated for DML localities using the ICE-D: DMC interface: http://dmc.ice-d.org by Perry Spector. Ice sheet models by Spector et al., 2018, Pollard & DeConto 2009, deBoer et al., 2014.
- The Function folder contains 
    - i) codes used for calculation of production rates retrieved (july 2021) from: https://bitbucket.org/cronusearth/cronus-calc/src/master/ ; Website: https://cronus.cosmogenicnuclides.rocks/2.1/ ; See Marrero et al., 2016 for details. 
    - ii) Muon production code by Balco et al., 2017
    - iii) export_fig package from https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig
    - iv) Sub-functions for visualisation purposes
- The model folder contains the sub-folder MDML used for storing model output. The MDML sub-folder reports stores reports generated by makereportJ_mn5.m

References
 - Balco, G. (2017) Production rate calculations for cosmic-ray-muon-produced 10Be and 26Al benchmarked against geological calibration data. Quaternary Geochronology, 39, 150-173.
 - De Boer, B., Lourens, L. J., & Van De Wal, R. S. (2014). Persistent 400,000-year variability of Antarctic ice volume and the carbon cycle is revealed throughout the Plio-Pleistocene. Nature communications, 5(1), 1-8.
 - Laskar, J., Fienga, A., Gastineau, M., & Manche, H. (2011). La2010: a new orbital solution for the long-term motion of the Earth. Astronomy & Astrophysics, 532, A89.
 - Lisiecki, L. E., & Raymo, M. E. (2005). A Pliocene‐Pleistocene stack of 57 globally distributed benthic δ18O records. Paleoceanography, 20(1).
 - Marrero, S. M., F. M. Phillips, B. Borchers, N. Lifton, R. Aumer & G. Balco (2016a) Cosmogenic nuclide systematics and the CRONUScalc program. Quaternary Geochronology, 31, 160-187.
 - Pollard, D., & DeConto, R. M. (2009). Modelling West Antarctic ice sheet growth and collapse through the past five million years. Nature, 458(7236), 329-332.
 - Spector, P., Stone, J., Pollard, D., Hillebrand, T., Lewis, C., & Gombiner, J. (2018). West Antarctic sites for subglacial drilling to test for past ice-sheet collapse. The Cryosphere, 12(8), 2741-2757.
 - Zachos, J., M. Pagani, L. Sloan, E. Thomas & K. Billups (2001) Trends, rhythms, and aberrations in global climate 65 Ma to present. Science, 292, 686-693.


/Jane Lund Andersen, Aarhus University, Sep. 2022
