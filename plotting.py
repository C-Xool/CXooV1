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

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.ticker as tick
from cartopy import feature as cfeature
import cmocean
import cmocean.plots
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

from specifications import ROMSSpecScalar, ROMSSpecVector

class ScalarFormat(matplotlib.ticker.ScalarFormatter):
    def __init__(self, fformat="%1.1f", offset=True, mathText=True):
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,
                                                        useMathText=mathText)

def get_cmap_limits(data):
    return ( np.floor(np.nanmin(data) ),
             np.ceil( np.nanmax(data) ) )
    
def get_ticks(limits, N):
    return np.round(np.linspace(*limits, N),2)

def get_data_source(spec, file):
    if isinstance(spec, ROMSSpecScalar):
        return file[spec.roms_name]
    elif isinstance(spec, ROMSSpecVector):
        return tuple( ( file[spec.roms_comp_names[0]], file[spec.roms_comp_names[1]] ) )
    else:
        raise Exception("incompatible spec.")

def plot_variables(var_list_to_plot,
                   loaded_file_to_plot,
                   discrete_colors,
                   homogenise_lims,
                   interval,
                   projected_coords,
                   projection,
                   var_angle,
                   scale_factor,
                   arrowdensity,
                   output_folder,
                   plots_folder
                  ):

    if not isinstance(var_list_to_plot, (list, tuple)):
        raise Exception("Variables have to be specified as list.")
    
    dates_to_plot=loaded_file_to_plot['time']
    
    for spec in var_list_to_plot:
        
        print(f"Plotting: {spec.roms_name}")
        
        geo_var = spec.roms_name
        
        containingfolder=os.path.join(output_folder,plots_folder,geo_var)

        os.makedirs(containingfolder,exist_ok=True)

        vars_to_plot = get_data_source(spec, loaded_file_to_plot)
        
        plot_single_variable(spec,
                             vars_to_plot,
                             dates_to_plot,
                             homogenise_lims,
                             interval,
                             projected_coords,
                             discrete_colors,
                             var_angle,
                             scale_factor,
                             arrowdensity,
                             containingfolder,
                             projection,
                             plots_folder
                        )

def data_to_plot( fun_geo_var, fun_time_input, fun_interval=None):
    if(fun_interval==None):
        slc = slice(None)
    elif(fun_interval>0):
        slc = slice(None, None, fun_interval)
    else:
        slc = slice(-fun_interval, -fun_interval+1)

    time_data = fun_time_input[slc]
    
    if isinstance(fun_geo_var, tuple):
        vector_data = tuple( arr.isel(time=slc) for arr in fun_geo_var )

        scalar_data = np.sqrt( sum( arr**2 for arr in vector_data) )
        scalar_data = scalar_data.transpose("time", "y",  "x")
    else:
        scalar_data = fun_geo_var.isel(time=slc)

        vector_data = None

    return time_data, scalar_data, vector_data


def plot_single_variable(spec,
                          var_to_plot_01,
                          dates_to_plot,
                          homogenise_limits,
                          interval,
                          projected_coord,
                          discrete_colors,
                          var_angle,
                          scale_factor,
                          arrowdensity,
                          containingfolder,
                          projection,
                          plots_folder
                        ):
    
    time_data, scalar_data, vector_data = data_to_plot(var_to_plot_01,
                                                       dates_to_plot,
                                                       interval)
  
    if(homogenise_limits):
        scalar_limits = get_cmap_limits(scalar_data)
    else:
        scalar_limits = None
                                                                   
    if ( not discrete_colors ):
        levelsgeovar = 256
    elif isinstance(discrete_colors, int) and discrete_colors > 1:
        levelsgeovar = discrete_colors + 1
    else:
        raise Exception(f"discrete values must be a positive integer larger than one or one of (0, False, None), got {discrete_colors}")

    if not discrete_colors:
        Nticks = 10
    else:
        Nticks = levelsgeovar

    if (vector_data is not None):        
        staffelberg=dict(eta_rho='y',
                         xi_rho='x')
                
        cosrot=np.cos(var_angle).rename(staffelberg)
        sinrot=np.sin(var_angle).rename(staffelberg)

        uv = xr.concat( [cosrot*vector_data[0] - sinrot*vector_data[1], 
                        sinrot*vector_data[0] + cosrot*vector_data[1] ], dim='uv')
        
    
        uv=uv.transpose('uv','time','y','x')

        
    if scalar_limits is not None:
        level_boundaries_geo_var = np.linspace(*scalar_limits, levelsgeovar)

        tick_values = get_ticks(scalar_limits, Nticks)
        

    for t, _date in enumerate(np.array(time_data)):
        
        geo_var = spec.roms_name
        
        file_name=str(geo_var)+'_date_'+str(_date)[0:19].replace('-','_').replace('T','_at_hr_').replace(':','_')+'.png'    

        if os.path.exists(containingfolder+'/'+file_name):
            print(str(file_name)+" already exits")
            continue

        fi = plt.figure(figsize=(24, 18), dpi=100)
        ax = plt.axes(projection=projection)
        ax.add_feature(cfeature.LAND.with_scale('50m'), linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=5, edgecolor='k')

        cmap = getattr(cmocean.cm, spec.plottype)

        if (scalar_limits is None):
            limits = get_cmap_limits(scalar_data[t])
            level_boundaries_geo_var = np.linspace(*limits, levelsgeovar)
            tick_values = get_ticks(limits, Nticks)

        s=ax.contourf(projected_coord[:,:,0], 
                      projected_coord[:,:,1],
                      scalar_data.isel(time=t), 
                      level_boundaries_geo_var, 
                      cmap=cmap)
             
        
        cb = plt.colorbar(s, ticks=tick_values,
                          boundaries=level_boundaries_geo_var,
                          values=(level_boundaries_geo_var[:-1] + level_boundaries_geo_var[1:]) / 2,
                          aspect=20, 
                          ax=plt.gca(),
                          alpha=0.75)
        
        cb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

        lname = spec.roms_long_name
        sname = spec.roms_name
        
        ax.set_title(lname+' ('+spec.plot_unit+')\n'+
                     str(_date)[0:19].replace('-',' ').replace('T',' at ')+' hr\n', fontsize=36)
        cb.set_label(sname+' ('+spec.plot_unit+')', size=28)
        
        if (vector_data is not None):
            

            ax.quiver(projected_coord[::arrowdensity,::arrowdensity,0],
                      projected_coord[::arrowdensity,::arrowdensity,1],
                      uv[0,t,::arrowdensity,::arrowdensity],
                      uv[1,t,::arrowdensity,::arrowdensity],
                      color='k',scale=scale_factor,
                      linewidths=1, width=0.0005, edgecolor='0', 
                      headwidth=7,headaxislength=7,headlength=15)
        
        cb.ax.tick_params(labelsize=20, length=0)
        cb.set_alpha(0.2)
        cb.outline.set_visible(False)
        grids = ax.gridlines(draw_labels=True, 
                             linewidth=1.0, 
                             color='gray', 
                             alpha=0.5, 
                             linestyle='--')
        
        grids.top_labels = True
        grids.bottom_labels = True
        grids.left_labels = True
        grids.right_labels = True
        grids.ylines = True
        grids.xlines = True  
        
        grids.xformatter = LongitudeFormatter()
        grids.yformatter = LatitudeFormatter()
        grids.xlabel_style = {'size': 15, 'color': 'gray'} 
        grids.ylabel_style = {'size': 15, 'color': 'gray'} 

        outputpath=os.path.join(containingfolder,file_name)
        print("Writing file: "+str(outputpath))
        fi.savefig(outputpath)
        
        plt.clf()
        plt.close()
        ax.clear()
