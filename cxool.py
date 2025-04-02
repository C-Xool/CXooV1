# -*- coding: utf-8 -*-
"""
CXool – V1: Oceanographic exploration, 
is a program to prepare the grid domain with the Atmospheric forcing,
to carry on with ocean modelling in ROMS. equations.

 -> This is a free software; you can redistribute it and/or
 -> modify it under the terms of the GNU General Public License
 -> as published by the Free Software Foundation; either version 3
 -> of the License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 -> but WITHOUT ANY WARRANTY; without even the implied warranty of
 -> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -> GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 -> along with this program.  If not, see <http://www.gnu.org/licenses/>.
 Author and main maintainer: Carlos Argáez
 -> Bibliography attached to the corresponding publication.
 Authors:
     Carlos Argáez, Simon Klüpfel, María Eugenia Allenda Aranda, Christian Mario Appendini
     To report bugs, questions, critics or just greetings, please use:
         cargaezg@iingen.unam.mx
 """

import os
import shutil
import argparse
import sys     
   
import numpy as np
import xarray as xr
from cartopy import crs as ccrs

from specifications import ERA5SpecPL, Rutgers as vardict
from cds_handler import CDSHandler, _known_pressure_levels
from interpolation import interpolate_to_grid
from plotting import plot_variables



class CXoolArgumentParser():
    def __init__(self, debug_args=[]):
        
        parser0 = argparse.ArgumentParser(
                    prog='C-Xool',
                    description='Oceanographic exploration - prepare the grid domain with the Atmospheric forcing, to carry on with ocean modelling in ROMS. equations',
                    epilog='To report bugs, questions, critics or just greetings, please use: cargaezg@iingen.unam.mx',
                    fromfile_prefix_chars='@',
                    add_help=False)
        
        argv0 = sys.argv[1:] + debug_args
        
        parser0.add_argument('-f','--config-file')
        
        cfa, argv1 = parser0.parse_known_args(argv0)
        if cfa.config_file is not None:
            argv = ['@'+cfa.config_file] + argv1
        else:
            argv = argv1

        del parser0
            
        parser = argparse.ArgumentParser(
                    prog='C-Xool',
                    description='Oceanographic exploration - prepare the grid domain with the Atmospheric forcing, to carry on with ocean modelling in ROMS. equations',
                    epilog='To report bugs, questions, critics or just greetings, please use: cargaezg@iingen.unam.mx',
                    fromfile_prefix_chars='@')

        parser.add_argument('-f','--input-file',
                            dest='config_file', 
                            help='Input file to read the instructions to run. It is generated as an external input file or created with the instruction <-o name-of-input-file>.')
        
        ##### THESE LINES ARE FOR GRID PARAMETERS
        parser.add_argument('-a', 
                            '--grid_name', 
                            required=True,
                            help='Name of the ROMS grid file.',
                            type=str)
        parser.add_argument('-b',
                            '--initialdate', 
                            required=True,
                            type=str,
                            help='Initial date that defines when the download and interpolation start. Format YYYY-MM-DD. The hour is automatically set as 00:00:00.') 
        parser.add_argument('-c', '--finaldate',  required=True,
                            help='Final date that defines when the download and interpolation end. Format YYYY-MM-DD. The hour is automatically set as 23:00:00.') 
        
        parser.add_argument('-e', '--final_interpolated_file', default = 'ERA5_interpolated_to_grid.nc',
                            type=str,
                            help="String, use \"\" (quotes) to define it. It is the name of the final interpolated NetCDF file with the atmospheric forcing.",
                            metavar='"<name_of_your_final_file.nc>"')
        ###### THESE LINES ARE FOR PLOTTING        
        parser.add_argument('-g', '--plot', 
                            dest='vars_to_plot',
                            nargs='*', 
                            default=[ ], 
                            choices=list(vardict),
                            type=str,
                            help="This command defines if C-Xool will plot variables from the final forcing file.",
                            metavar='<var1 var2...>')   #Use the world wind for variable winds
        parser.add_argument('-gg', '--vars_to_interp', nargs='+', 
                            default=[ 'msl', 'wind', 't2m','tcc'], 
                            type=str,
                            help="Defines the variables to download and to interpolate. The default variables are: msl, wind, t2m, tcc.",
                            metavar='<var1 var2...>')   #Use the world wind for variable winds
        parser.add_argument('-g2', '--projection', default="mercator", 
                            choices=["stereographic", "mercator"], 
                            type=str,
                            help='Which map projection is used for plotting the interpolated data out of the final netcdf file.' )
        
        parser.add_argument('-i', '--interval', default=6, type=int, 
                            help='Defines the time interval, in hours, to download the data: Every <interval> hours counting from 00:00 on the given date.', 
                            metavar='<Number of intervals between hours>')
        parser.add_argument('-i2', '--plot_interval', default=10, type=int, 
                            help='Defines the time interval, in hours, to plot the interpolated of the final netcdf file: Every <interval> hours counting from 00:00 on the given date.', 
                            metavar='<Number of intervals between hours>')
        parser.add_argument('-j', '--scale_factor',
                            help='Defines the scale, in Python syntax, used by the quiver plot to control the arrow size. It does the scaling length of the arrow inversely. Default 200.',
                            metavar='<Value>',
                            default=200, 
                            type=float)
        parser.add_argument('-k', '--arrowdensity', 
                            default=1, 
                            type=int,
                            help='Used by C-Xool to defined the density of the arrows to be plotted. it does it inversely: The highest the value input the less arrows to be seen.',
                            metavar='<Value>')
        parser.add_argument('-l', '--discrete_colors', 
                            default=10, 
                            type=int,
                            help='For analysis: The discretisation of colours allows to clearly see the fronts in the atmospheric variables. Chose not to use would be the colours change is soft.',
                            metavar='<Value>')
        parser.add_argument('-m','--homogenise_limits', 
                            action=argparse.BooleanOptionalAction, 
                            default=True,
                            help='For analysis: Homoegenising the limits permits that all the plotted variables have their colour bar from <-extreme_value> to <extreme_value>, permitting a better analysis.`',
                            metavar='<var1 var2...>')
        parser.add_argument('-o', '--out-config', 
                            default = None,
                            help="Generates an input file to be read. It is recommended to create it.",
                            metavar="<input_file_name>",
                            type=str)
        parser.add_argument('-n', '--output_folder', default = '.',
                            type=str,
                            help="String, use \"\" (quotes) to define it. It is the name of the folfer where the output will be saved.",
                            metavar='"Folder_name"')
        parser.add_argument('-p', '--data_storage', default = None,
                            type=str,
                            help="String, use \"\" (quotes) to define it. It is the name of the folfer where the ERA5 raw data will be saved.",
                            metavar='"/path/to/data/storage"')
        parser.add_argument('-q', '--plots_folder', default = "plots",
                            type=str,
                            help="String, use \"\" (quotes) to define it. It is the name of the subfolfer where the plots will be saved.",
                            metavar='<"name_of_folder">')
        parser.add_argument('-r', '--data_subfolder', default = "merged_data",
                            type=str,
                            help="String, use \"\" (quotes) to define it. It is the name of the subfolfer where the raw merge data will be saved.",
                            metavar='<"name_of_folder">')
        parser.add_argument('-s', '--memory_chunks', 
                            default=None, 
                            type=int,
                            help='---.',
                            metavar='<Value>')
        


        
        


        
        self.parser = parser
        self.args = self.parser.parse_args(argv, )
        self.argv = argv
        
        if self.args.out_config is not None and cfa.config_file == self.args.out_config:
            raise Exception('Input config file would be overwritten. Aborting.')
            
        if not set(self.args.vars_to_plot).issubset(set(self.args.vars_to_interp)):
            raise Exception("Your lists for both interpolating and plotting do not match.")
        
        if self.args.out_config is not None:
            dels = []
            for iv, v in enumerate(argv0):
                if v in ('-o','--out_config'):
                    dels.extend([iv, iv+1])
                dels = sorted(set(dels),reverse=True)
                for d in dels:
                    del argv0[d]

            with open(self.args.out_config,'wt') as of:
                of.writelines('\n'.join(argv0))
                
                


def padl(x):
    return np.floor(x*4.0)/4.0-0.5
def padr(x):
    return np.ceil(x*4.0)/4.0+0.5

class CXool():
    def __init__(self, grid, 
                 ini_date, 
                 fin_date, 
                 interval, 
                 variable_list, 
                 data_storage,
                 output_folder,
                 plots_folder,
                 data_subfolder,
                 memory_chunks,
                 out_config
                 ):
        self.data_storage = data_storage
        self.output_folder=output_folder
        self.plots_folder=plots_folder
        self.data_subfolder=data_subfolder
        self.out_config=out_config
        self.interval = interval
        self.grid = grid.copy()
        lat = self.grid.lat_rho.values 
        lon = self.grid.lon_rho.values
        self.ini_date=ini_date
        self.fin_date=fin_date

        self.variable_list=variable_list
        
        self.lim_lat = ( lat.min(), lat.max() )
        self.lim_lon = ( lon.min(), lon.max() )
        
        self.lim_lat_pad = ( padl(self.lim_lat[0]), padr(self.lim_lat[1]) )
        self.lim_lon_pad = ( padl(self.lim_lon[0]), padr(self.lim_lon[1]) )

        self.lim_NWSE = ( lat.max(), lon.min(), lat.min(), lon.max() )
        padfuns = ( padr, padl, padl, padr )
        self.lim_NWSE_pad = tuple( ( pf(x).item() for pf, x in zip(padfuns, self.lim_NWSE) ) )
        
        self._ts_paths = None
        self.memory_chunks=memory_chunks
    def move_out_config(self):
        if self.out_config is not None:
            folder_path = self.output_folder
            os.makedirs(folder_path, exist_ok=True)
            shutil.move(self.out_config, folder_path)


        
    def prepare_timeseries(self, force_splice = False, force_download = False):
        
        if not force_splice and self._ts_paths is not None:
            return
        self._ts_paths = {}


        CDSH = CDSHandler(self.data_storage)

        vlist = sum( ( v.era5spec for v in self.variable_list.values() ), start=[] )
        
        for variable in vlist:
            varlab = variable.name
            print(f'Processing raw data for variable: {varlab}')
            fname=f'{self.output_folder}/{self.data_subfolder}/{varlab}_'+str(self.ini_date)+'_'+str(self.fin_date)+'.nc'
            
            self._ts_paths[varlab] = fname

            if isinstance(variable, ERA5SpecPL):
                pressure_level = variable.pressure_level
            else:
                pressure_level = None
                
            tsds = CDSH.get_timeseries(variable.long_name,
                                       self.ini_date, 
                                       self.fin_date, 
                                       self.lim_NWSE_pad, 
                                       self.interval,
                                       pressure_level
                                       )
            os.makedirs(os.path.dirname(fname),exist_ok=True)
        
            if tsds[varlab].isnull().sum().compute().item() > 0:
                raise Exception('Nans have been found in ERA5 data. Execution stopped.')

            
            tsds.to_netcdf(fname) # Export netcdf file and save it.
            tsds.close()

    def interpolate(self, filename):
        interpolate_to_grid(self.grid, 
                            filename, 
                            self.variable_list, 
                            self._ts_paths,
                            self.memory_chunks) 



    def plotit(self,
               file_to_plot, 
               gridfile, 
               varlist,
               projection='mercator',
               interval=10,
               scale_factor=500,
               arrowdensity=15,
               discrete_colors=10,
               homogenise_limits=False,
               output_folder='.',
               plots_folder=None
               ):
        
        loaded_file_to_plot = xr.open_dataset(file_to_plot)
    
        with xr.open_dataset(gridfile) as new:
            angle = new.angle
            cen_lat = float(((new.lat_rho.max()+new.lat_rho.min())/2))
            cen_lon = float(((new.lon_rho.max()+new.lon_rho.min())/2))
    
        wgs84 = ccrs.CRS('WGS84', globe=None)
        
        if (projection=='stereographic'):
            projection = ccrs.Stereographic(central_latitude=cen_lat, 
                                            central_longitude=cen_lon, 
                                            false_easting=0.0, 
                                            false_northing=0.0, 
                                            true_scale_latitude=None, 
                                            globe=None
                                            )
            
        elif (projection=='mercator'):
            projection = ccrs.Mercator(central_longitude=cen_lon,
                                  false_easting=0.0, false_northing=0.0,
                                  latitude_true_scale=None, globe=None)
        else:
            raise Exception(f'unknown projection: "{projection}"')
        
        x = loaded_file_to_plot['lon']
        y = loaded_file_to_plot['lat']
        xx = projection.transform_points(wgs84, 
                                         np.array(x).flatten(), 
                                         np.array(y).flatten()
                                         )
    
        projected_grid = xx.reshape( x.shape + (3,) )
    
        plot_variables(varlist,
                       loaded_file_to_plot,
                       discrete_colors,
                       homogenise_limits,
                       interval,
                       projected_grid,
                       projection,
                       angle,
                       scale_factor,
                       arrowdensity,
                       output_folder,
                       plots_folder)
    

dbgargs = []

if __name__ == "__main__":
        
    try:
        CXA = CXoolArgumentParser(dbgargs) 
    except Exception as e:
        print('Error: ', str(e))
        sys.exit(1)
    
    instructions = CXA.args

    ps = { a.split('/')[0]:(int(a.split('/')[1]) if '/' in a else None) for a in instructions.vars_to_interp }
        
    for v in ps.keys():
        if v not in vardict:
            raise KeyError(f"Wrong variable name: The {v} is not a valid variable.")

    for p in ps.values():
        if p is not None and p not in _known_pressure_levels:
            raise ValueError(f"Wrong level pressure value: the level {p} is not a valid value.")
    
    vd = { a:vardict[a] for a in ps }
    for v, p in ps.items():
        if p is not None:
            vd[v].era5spec[0].pressure_level=p

    
    grid = xr.open_dataset(instructions.grid_name)
    
    CX = CXool(grid, 
               instructions.initialdate, 
               instructions.finaldate, 
               instructions.interval, 
               vd, 
               instructions.data_storage,
               instructions.output_folder,
               instructions.plots_folder,
               instructions.data_subfolder,
               instructions.memory_chunks,
               instructions.out_config)
    
    
    CX.prepare_timeseries()
    ncf_filename=os.path.join(CX.output_folder,instructions.final_interpolated_file)

    CX.interpolate(ncf_filename)
    
    CX.move_out_config()
      
    if len(instructions.vars_to_plot):

        varlist = [ vd[v] for v in instructions.vars_to_plot ]
        CX.plotit(ncf_filename,
                  instructions.grid_name,
                  varlist,
                  instructions.projection,
                  instructions.plot_interval,
                  instructions.scale_factor,
                  instructions.arrowdensity,
                  instructions.discrete_colors,
                  instructions.homogenise_limits,
                  instructions.output_folder,
                  instructions.plots_folder
                  )
            
    print("Task accomplished, bye.")


