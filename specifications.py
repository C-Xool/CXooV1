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
from dataclasses import dataclass

@dataclass
class ERA5SpecSL:
    """Class for keeping track of the variable specs and scalings."""
    name: str
    long_name: str
    unit: str

@dataclass
class ERA5SpecPL(ERA5SpecSL):
    pressure_level: int


    

@dataclass
class ROMSSpecScalar():
    era5spec: list
    #
    roms_name: str
    roms_long_name: str
    roms_unit: str
    plot_unit: str
    #transform
    scale: float = 1.0
    shift: float = 0.0
    plottype: str = None
    
    def __post_init__(self):
        assert isinstance(self.era5spec, ERA5SpecSL)
        self.era5spec = [ self.era5spec ]

@dataclass
class ROMSSpecVector():
    era5spec: list
    #
    roms_name: str
    roms_long_name: str
    roms_unit: str
    #transform
    roms_comp_names: list
    roms_comp_long_names: list
    plot_unit: str
    scale: float = 1.0
    shift: float = 0.0
    plottype: str = None
    

    @property
    def name(self):
        raise Exception('nonono')

    def __post_init__(self):
        assert len(self.era5spec) == 2


Rutgers ={'msl' : ROMSSpecScalar( 
                   ERA5SpecSL('msl', 'mean_sea_level_pressure', 'Pa'),
                   'Pair', 'Mean sea level pressure', 'hPa', 'hPa',
                   scale=1.0e-2, plottype='diff' ),
           't2m' : ROMSSpecScalar( 
                   ERA5SpecSL('t2m', '2m_temperature', 'K'),
                   'Tair', '2 metre temperature', 'C', '°C',
                   shift=-273.15, plottype='thermal'),
           #Qair or total_cloud_cover are not available in ERA5, only available in reanalysis-era5-pressure-levels.
           'tcc' : ROMSSpecScalar( 
                   ERA5SpecSL('tcc', 'total_cloud_cover', '(0-1)'),
                   'cloud', 'Total cloud cover', '(0-1)', '(0-1)',
                   plottype='balance'),
           'avg_snswrf' : ROMSSpecScalar( 
                      ERA5SpecSL('avg_snswrf', 'mean_surface_net_short_wave_radiation_flux', 'W m-2'),
                               'swrad', 'Mean surface net short-wave radiation flux', 'W m-2', 'W m$^2$',
                               plottype='solar'),
           'avg_sdlwrf': ROMSSpecScalar( 
                       ERA5SpecSL('avg_sdlwrf', 'mean_surface_downward_long_wave_radiation_flux', 'W m-2'),
                               'lwrad_down', 'Mean surface downward long-wave radiation flux', 'W m-2', 'W m$^2$',
                               plottype='solar'),
           'avg_snlwrf' : ROMSSpecScalar( 
                       ERA5SpecSL('avg_snlwrf', 'mean_surface_net_long_wave_radiation_flux', 'W m-2'),
                               'lwrad', 'Mean surface net long-wave radiation flux','W m-2', 'W m$^2$',
                               plottype='solar'),
           'avg_slhtf' :   ROMSSpecScalar( 
                       ERA5SpecSL('avg_slhtf','mean_surface_latent_heat_flux', 'W m-2'),
                               'latent', 'Mean surface latent heat flux', 'W m-2', 'W m$^2$',
                               plottype='solar'),
           'avg_ishf' :   ROMSSpecScalar( 
                       ERA5SpecSL('avg_ishf', 'mean_surface_sensible_heat_flux', 'W m-2'),
                               'sensible', 'Mean surface sensible heat flux','W m-2', 'W m$^2$',
                               plottype='solar'),
           'tp'  :     ROMSSpecScalar( 
                       ERA5SpecSL('tp', 'total_precipitation', 'm'),
                               'rain', 'Total precipitation', 'kg m-2 s-1', 'kg m$^{-2}$ s$^{-1}$',
                               scale=10.0/36.0, plottype='rain'),
           'e'   :     ROMSSpecScalar( 
                       ERA5SpecSL('e', 'evaporation', 'm of water equivalent'),
                               'evap', 'Evaporation', 'kg m-3','kg m$^{-3}$',
                               scale=10.0/36.0,
                               plottype='rain'),
           'cc' : ROMSSpecScalar( 
                       ERA5SpecPL('cc', 'fraction_of_cloud_cover', '0-1', pressure_level=300),
                               'fcloud', 'cloud fraction', '0-1', '0-1',
                               plottype='balance'),
           'q'   : ROMSSpecScalar( 
                       ERA5SpecPL('q', 'specific_humidity', 'kg kg-1', pressure_level=1000),
                               'Qair', 'surface specific humidity', 'kg kg-1', 'kg kg${-1}$',
                               scale=10.0/36.0,
                               plottype='rain'),
           'wind' :     ROMSSpecVector( 
                           [ ERA5SpecSL('u10', '10m_u_component_of_wind', 'ms-1'),
                             ERA5SpecSL('v10', '10m_v_component_of_wind', 'ms-1'),
                           ] ,

                        roms_name='Wind',
                        roms_long_name='10 metre wind',
                        roms_unit='ms-1',
                        plot_unit='ms$^{-1}$',
                        roms_comp_names = ['Uwind', 'Vwind'],
                        roms_comp_long_names = ['10 metre U wind component',
                                                '10 metre V wind component'],
                        
                        plottype='speed'

                        )
           }

