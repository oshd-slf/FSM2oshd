# FSM2oshd

FSM2oshd (Mott et al., under review) solves the complete mass and energy balance of the snowpack for open and forested areas and is originally based on the Flexible Snow Model FSM2 [(Essery, 2015)](#Essery2015). Specifically designed for snow-hydrological forecasting in Switzerland, FSM2oshd is developed by the operational snow-hydrological service at the WSL Institute for Snow and Avalanche Research SLF [operational snow-hydrological service](https://www.slf.ch/en/snow/snow-as-a-water-resource/snow-hydrological-forecasting.html). 


FSM2oshd includes novel model components that account for forest processes [(Mazzotti et al., 2021)](#Mazzotti2021) and snow cover fraction [(Helbig et al., 2021)](#Helbig2021) computations. For larger-scale applications sub-grid parameterizations for radiation [(Helbig and Löwe, 2014)](#Helbig2014) and wind [(Helbig et al., 2017)](#Helbig2017) have been implemented. Additionally, FSM2oshd features updated snow pack process parameterizations, including snow compaction, new snow density, decay of snow albedo, and snow hydraulics. 

In addition to its point-based capabilities, FSM2oshd can also be run on a regular grid. When run on a grid, the model divides each grid cell into separate tiles, accounting for forest-covered, open, and glacierized fractions. FSM2oshd, in its grid mode, includes several enhancements for the forest tiles [Mazzotti et al., (2020)](#Mazzotti2020). These improvements include the incorporation of local canopy cover fraction, hemispherical sky view fraction, wind attenuation, and spatial patterns of shortwave and longwave radiation transfer through the canopy [(Mazotti et al., 2023)](#Mazzotti2023). Furthermore, the model now accounts for preferential deposition of snow in canopy gaps, resulting in more accurate predictions. These enhancements ensure that FSM2oshd provides highly detailed and precise predictions for snow-hydrological processes in complex forested terrain.

The FSM2oshd model can be utilized for seasonal runs to analyze and predict spatio-temporal dynamics of snow processes. Output variables include snow water equivalent (SWE), snow depth (HS), snow melt runoff (SMR), and snow state, based on liquid water content. The model provides predictions at both daily and hourly intervals, allowing for a comprehensive understanding of the temporal evolution of snow-related variables. 

## Building and running the model

The main computation routines of FSM2oshd are located in the **./src/core** folder. For demonstration purposes, the model
can be run for a point simulation using input data stored in a text file (**./data**). In addition to the meteorological input data used in the original FSM2 version [(Essery, 2015)](#Essery2015), FSM2oshd requires input of 24h snowfall sums and canopy transmissivity for direct shortwave radiation (which should be approximated by the constant value of sky view fraction if such a time series is not available). The code for handling model setup and input/output operations related to the text file driving option is found in the **./src/txt** folder. 

We are currently working on creating example files that can be used to run FSM2oshd in its grid mode, and we anticipate making them available in the near future.

### Windows compilation
To compile the model for windows, use the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler and execute the compile script in the root folder:

`./compil_txt.bat`

Run the executable with the command:

`./FSM_TXT.exe nlst_txt.nam`

Model output is stored in a text file (see **./data**).

### Linux compilation

```
#-----------------------
#  compile using gfortran
#-----------------------

bash compile_oshd_txt.sh

#-----------------------
#  execute example
#-----------------------

# open mode
./FSM_OSHD nlst_txt.nam

# forest mode
./FSM_OSHD nlst_txt_for.nam
```

### Python wrapper

See the Python script `runner.py` for performing simulations and displaying the results for different model configurations (albedo, canopy, conductivity, snow density, hydraulics). The `requirements.txt` file contains information about the virtual environment used for running this code. The file nlst_txt_forest.nam contains the settings adapted to simulations at forest points.


## References


<a name="Essery2015"></a> Essery (2015). A Factorial Snowpack Model (FSM 1.0). *Geoscientific Model Development*, **8**, 3867-3876, [doi:10.5194/gmd-8-3867-2015](http://www.geosci-model-dev.net/8/3867/2015/)

<a name="Helbig2014"></a> Helbig, N., & Löwe, H. (2014). Parameterization of the spatially averaged sky view factor in complex topography. Journal of Geophysical Research: Atmospheres, 119(8), 4616-4625.

<a name="Helbig2017"></a> Helbig, N., Mott, R., Van Herwijnen, A., Winstral, A., & Jonas, T. (2017). Parameterizing surface wind speed over complex topography. Journal of Geophysical Research: Atmospheres, 122(2), 651-667.

<a name="Helbig2021"></a> Helbig, N., Bühler, Y., Eberhard, L., Deschamps-Berger, C., Gascoin, S., Dumont, M., Revuelto, J., Deems, J. S., and Jonas, T. (2021a). Fractional snow-covered area: scale-independent peak of winter parameterization, The Cryosphere, 15, 615–632, https://doi.org/10.5194/tc-15-615-2021.

<a name="Mazzotti2020"></a> Mazzotti, G., Essery, R., Moeser, C.D. and Jonas, T. (2020). Resolving small‐scale forest snow patterns using an energy‐balance snow model with a 1‐layer canopy. Water Resources Research, 56(1), e2019WR026129. [doi: 10.1029/2019WR026129 ](https://)

<a name="Mazzotti2021"></a> Mazzotti, G., Webster, C., Essery, R. and Jonas, T. (2021) Increasing the physical representation of forest snow processes in coarse-resolution models: lessons learnt from upscaling hyper-resolution simulations. Water Resources Research 57(5), e2020WR029064. [doi: 10.1029/2020WR029064.](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020WR029064)

<a name="Mazzotti2023"></a>Mazzotti, G., Webster, C., Quéno, L., Cluzet, B., and Jonas, T. (2023). Canopy structure, topography and weather are equally important drivers of small-scale snow cover dynamics in sub-alpine forests, Hydrol. Earth Syst. Sci. Discuss. [preprint], https://doi.org/10.5194/hess-2022-273, in press. 




