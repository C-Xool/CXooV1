<div align="right">
This file was written at LIPC, IINGEN, UNAM, <br>
in Sisal, Yucatan, Mexico, <br>
on the 9th of March 2025.
</div>

# C-Xool Software Documentation

**This Software is called C-Xool and its distribution and modifications**  
**are in accordance with the GNU GENERAL PUBLIC LICENSE.**
  
This file is offered as-is, without any warranty.


# Table of Contents

1. [Introduction](#introduction)
2. [Developers](#developers)
3. [Installing C-Xool](#installing-c-xool)
4. [Architecture](#architecture)
5. [Execution](#execution)
    - [Basic Execution](#basic-execution)
    - [Commands](#commands)
    - [Example of execution](#example-of-execution)
    - [Help and Further Capabilities](#help-and-further-capabilities)
6. [Support](#support)
7. [Known Problems](#known-problems)
8. [Platforms](#platforms)
    - [Note for Mac OSX Users](#note-for-mac-osx-users)
    - [Note for Windows Users](#note-for-windows-users)



## Introduction


C-Xool is a toolbox designed to provide a powerful tool for computational analysis and data handling for ocean modelling within a framework of Earth Sciences. This is the first-ever implementation, specifically designed to build atmospheric forcing. It offers a more efficient and flexible solution for processing data for ocean modelling and related applications.

C-Xool was designed and developed at The Coastal Engineering and Processes Laboratory (LIPC, from its Spanish acronym *Laboratorio de Ingeniería y Procesos Costeros*), at the Engineering Institute of the National Autonomous University of Mexico (UNAM, *Universidad Nacional Autónoma de México*), in an office that overlooks the Gulf of Mexico and its beautiful landscape.

C-Xool is freely available for download. You are welcome to clone, modify, and share it with your friends and colleagues, as long as you comply with the terms of the GNU GENERAL PUBLIC LICENSE.

We encourage contributions and improvements, so if you have any suggestions or would like to share your work, feel free to reach out.

Please do not hesitate to send your comments, feedback, or any questions to cargaezg@iingen.unam.mx. We appreciate your input!


---

## Developers

The developers of C-Xool are:
  
  - Carlos Argáez García, LIPC, IINGEN, UNAM. 
  - Simon Klüpfel, Axelyf, Hafnarfjörður, Iceland
  - María Eugenia Allende Arandia, LIPC, IINGEN, UNAM. 
  - Christian Mario Appendini Albrechtsen, LIPC, IINGEN, UNAM. 


---

## Installing C-Xool


No installation needed. This code runs in Python and its execution should be straightforward. Please follow the general instructions:

1. Install a Python environment where to run C-Xool.
2. Install the following Python libraries:

    - cartopy  
    - cdsapi  
    - cmocean  
    - copernicusmarine  
    - dask  
    - matplotlib  
    - netcdf4  
    - numpy  
    - xarray  
    - argparse  


Once you have installed all the required libraries, enjoy C-Xool!

### Installation of the Mamba environment

````
mamba create -n cxool matplotlib scipy numpy netcdf4 dask xarray cdsapi cartopy cmocean
````



---


## Architecture:

1. C-Xool contains four scripts:

    - cxool.py<br>It contains the entry point of the executable and it is the neurological centre of the execution. 
    - cds_handler.py<br>
    This file handles the API request to Mercartor to download the data.	
    - interpolation.py<br>
    This file contains the functions to interpolate the atmospheric data onto the grid.
    - plotting.py<br>
    This file handles the plotting routines.
    - specifications.py<br>
    It contains variables dictionaries and variables classes.
---



## Execution
### Basic Execution
C-Xool runs on a Linux (or Python) terminal. It is given commands that allow controlling the execution parameters:

### Commands
Commands are introduced through the terminal. They are divided into mandatory and needed commands, and optional commands.

 **Needed commands:**  

**grid_name**  
String command. It is the name the user has given to their model grid, where the smoothed bathymetry to the model is located.  

**initialdate**  
Numeric command. It is the initial date that C-Xool will take to download the data to interpolate into the grid. The initial time is taken to be 00:00:00 hrs on the initial date.  

**finaldate**  
Numeric command. It is the final date that C-Xool will take to download the data to interpolate into the grid. The final time is taken to be 23:00:00 hrs on the final date.  


 **Optional commands:**  

**interval**  
Numeric command. It defines the time step to be used to download the data. Starting from the initial date at 00:00:00 hrs, the data is downloaded every *interval* hours until the final date is reached.  

**vars_to_interp**  
String commands. These are the names, in ERA5 syntax, of the atmospheric variables one needs to download. However, one variable is referred to differently: *wind*. This decision was made by the developers to ensure that both components, u and v, are properly downloaded and processed up to the completion of the interpolated file.  

**final_interpolated_file**  
String command. It is the name of the final atmospheric forcing file to be fed to the model.  

**plot_interval**  
String commands. It contains the variable names to be plotted. The user must ensure these names follow the ERA5 naming convention, except for wind.  

**projection**  
String command. In the current version of C-Xool, it can be either stereographic or mercator.  

**scale_factor**  
Numeric command. It controls the scale factor of the arrows in the quiver plot (vectorial variables).  

**arrowdensity**  
Numeric command. It controls the density of the arrows in the quiver plot (vectorial variables). A lower value results in higher arrow density.  

**discrete_colors**  
Numeric command. It generates a contour plot with clear colour distinctions to highlight different fronts. However, if set to None, it generates a plot with continuous colours.  

**homogenise_limits**  
Boolean command. If set to True, it generates a set of plots for a given variable with the same limits across the set. This is useful when comparing the evolution of the input variable over a period of time.  

**data_storage**  
String command. It is the name of the folder where the ERA5 raw data will be saved.  

**data_subfolder**  
String command. It is the name of the subfolder where the raw merged data will be saved.  

**plots_folder**  
String command. It is the name of the subfolder where the plots will be saved.  

**output_folder**  
String command. It is the name of the folder where the output will be saved.  

**final_interpolated_file**  
String command. It is the name of the final interpolated NetCDF file containing the atmospheric forcing.  

**memory_chunks**
Numeric command. When working on limited RAM-memory computers, it allows to divide the tasks into smaller chunks. 
It may reduce the speed of execution, but in turn, it will allow the process to complete. It reads an integer value, {e.g.} 100.

---

### Example of execution

#### To execute, run the following command:

```
python cxool.py --grid_name model_grid.nc --initialdate 1983-10-25 --finaldate 1983-10-27 
--vars_to_interp wind t2m msl tcc --plot msl wind t2m tcc --plot_interval 3
--final_interpolated_file="ERA5_interpolated_to_grid.nc" 
--output_folder "SoftwareXExample" -o "input"
```

---

### Help and further capabilities
#### To execute for help and further capabilities, run the following command:

```
python cxool.py --help
```

This will give you a list of functionalities that you can use out of C-Xool.

---
## Support
**In case of trouble or doubts, please send a description of the problem to
<cargaezg@iingen.unam.mx>.
You are welcome to contact the developer in English, Spanish, Italian or Yucatec Mayan.**

---

## Known problems

If you do not properly install the [Copernicus Python libraries](https://help.marine.copernicus.eu/en/articles/4854800-how-to-open-and-visualize-copernicus-marine-data-using-python) and set up your account on the [Copernicus](https://www.copernicus.eu/en) webpage, C-Xool may fail to execute.

If you do not properly install all the main needed Python libraries, as described above, C-Xool will not execute.

C-Xool requires as input:

    - A grid.
    - An initial date.
    - A final date.

#### If you do not provide these, C-Xool will fail to execute. All other commands are optional. 
---

## Platforms

### Note for Mac OSX users:

This code has been developed in Linux, which like OSX is a Unix based operating systems. The code has been tested on OSX and works. In case of any doubt, or error in its execution or compilation, please do contact the developer at <cargaezg@iingen.unam.mx>.


### Note for Windows users:

This code has been tested, used and partially developed using Windows, within a Python environment.

---
