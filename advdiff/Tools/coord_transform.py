# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 16:09:17 2022

This module contains functions for transforming xy UTM coordinates to lonlat coordinates and vice versa.

@author: Ketil
"""
import logging
import configparser
import numpy  as np
import pyproj as pp
import xarray as xr


def add_UTM_coords(dataset, store_AOI=True, risk=False, coords={'x'   : 'x',
                                                                'y'   : 'y',
                                                                'lon' : 'lon',
                                                                'lat' : 'lat'}):
    """
    Add UTM x-y coordinates to the dataset if they are not already there.
    Converts from longitude latitude coordinates using the coord_converter class.

    Parameters
    ----------
    dataset   : xarray dataset
                An xarray dataset containing latitude longitude coordinates to be converted into x-y UTM coordinates.
    store_AOI : bool
                Whether or not to store an area of interest file for fetching UTM zone for future coordinate transforms.
    risk      : bool
                If dataset we are adding UTM coordinates to are of standard format or riskmap format. 
    coords    : dict
                Dictionary containing x-y and lon-lat coordinate names.

    Returns
    -------
    xarray dataset
            An xarray dataset containing x-y UTM coordinates in coordinate.
    """
    # If x-y UTM coordinates not in dataset, convert from lon-lat to x-y.
    if not all(var in list(dataset.keys()) + list(dataset.coords) for var in [coords['x'],coords['y']]): 
        lonmin = dataset[coords['lon']].min().data
        lonmax = dataset[coords['lon']].max().data
        latmin = dataset[coords['lat']].min().data
        latmax = dataset[coords['lat']].max().data
        cc = coord_converter((lonmin, lonmax, latmin, latmax), store_AOI=store_AOI)

        if not risk: # All other files can easily be converted from lonlat to UTM
            dataset = cc.assign_lonlat_to_xy(dataset)
        else: # Riskmaps need some workarounds
            x, y = cc.lonlat_to_xy(dataset[coords['lon']].data.flatten(), dataset[coords['lat']].data.flatten())
            dataset = dataset.rename_dims({'lon':'x', 
                                           'lat':'y'})
            dataset = dataset.rename_vars({coords['lon']:coords['x'], 
                                           coords['lat']:coords['y']})
            dataset = dataset.assign({coords['x']: (['y','x'], x.reshape(dataset[coords['x']].data.shape)),
                                      coords['y']: (['y','x'], y.reshape(dataset[coords['y']].data.shape))})
    return dataset


class coord_converter:

    def __init__(self, lonlatbounds=None, store_AOI=True):
        """
        A coordinate converter class. This object can convert between x,y UTM coordinates and lon,lat coordinates on the Earth.
        For accurate results, longitude/latitude bounds need to be given. They are used to determine the relevant UTM zone.
        This is because UTM coordinates are NOT unique, so going from UTM to lon,lat is not defined without a UTM zone, however
        a UTM zone can be computed if you give a rough estimate for the lon,lat coordinates you are looking for with the lonlatbounds variable.
        If bounds are not given, the object will look for a UTM_AreaOfInterest.nc netcdf file which can give sufficient bounds.
        If you give lonlatbounds, a UTM_AreaOfInterest.nc file will be created such that you can call the converter again without having to give the bounds again.
        If that file has not been generated in the current run, the object will set the bounds to lon = [-180,180], lat  =[-90,90].
        However, this may not give accurate results as the UTM zones wil be wrong.

        Parameters
        ----------
        lonlatbounds : array ([lonmin, lonmax, latmin, latmax])
                       An array of the form [lonmin, lonmax, latmin, latmax].
        store_AOI    : bool
                       Whether or not to store the UTM_AreaOfInterest netcdf file whenever lonlatbounds are given.
        """
        config=configparser.ConfigParser(allow_no_value=True)
        config.read('External-Indata/AdvDiff.ini')
        self.datum  = config['coordinates']['datum']

        in_path  = config['paths']['indata_path'] + config['paths']['coord_path']
        filename = 'UTM_AreaOfInterest.nc'
        
        # We need to know the approximate bounds for the latitude and longitude coordinates to get the correct UTM zone. 
        # On the init of the first instance, the lonlat bounds will be given and they are stored in UTM_AreaOfInterest.nc. 
        # For other inits lonlat bounds may not be given and they will instead be loaded from UTM_AreaOfInterest.nc
        # As such, for every run, coordinate transforms will be consistent.
        if lonlatbounds is not None and store_AOI: # If bounds are given, store area of interest to file.
            self.bounds = lonlatbounds
            tmp = xr.Dataset()
            tmp.coords['lon']=np.array([self.bounds[0], self.bounds[1]])
            tmp.coords['lat']=np.array([self.bounds[2], self.bounds[3]])
            tmp.to_netcdf(in_path+filename) 
            tmp.close()                     
        else: # Else, try to load area of interest from file.
            try: # Try to load UTM_AreaOfInterest.nc
                tmp = xr.open_dataset(in_path+filename).load()
                tmp.close()
                lonmin = tmp.lon.min().data
                lonmax = tmp.lon.max().data
                latmin = tmp.lat.min().data
                latmax = tmp.lat.max().data
            
            except: # If UTM_AreaOfInterest.nc is not found, default. (This will give wrong coordinate transforms)
                logging.critical('No area of interest file found to convert from UTM to longitude latitude coordinates without UTM zone.')
                logging.critical('Will assume an area of interest of lon=[-180,180] and lat=[-90,90].')
                logging.critical('Wrong coordinate transforms will occur...')
                lonmin = -180.
                lonmax =  180.
                latmin = -90.
                latmax =  90.
                tmp = xr.Dataset()
                tmp.coords['lon']=np.array([lonmin, lonmax])
                tmp.coords['lat']=np.array([latmin, latmax])
                tmp.to_netcdf(in_path+filename) 
                tmp.close() 

            self.bounds = (lonmin, lonmax, latmin, latmax)

        lonmin, lonmax, latmin, latmax = self.bounds[:]
        self.area_of_interest = pp.aoi.AreaOfInterest(west_lon_degree  = lonmin,
                                                      south_lat_degree = latmin,
                                                      east_lon_degree  = lonmax,
                                                      north_lat_degree = latmax)
        epsg_code = pp.database.query_utm_crs_info(datum_name=self.datum, area_of_interest=self.area_of_interest)[0].code
        utm_crs   = pp.CRS.from_epsg(epsg_code)

        self.lonlat_xy_proj = pp.Transformer.from_crs(utm_crs.geodetic_crs, utm_crs, always_xy=True) # Transformer lonlat -> UTM
        self.xy_lonlat_proj = pp.Transformer.from_crs(utm_crs, utm_crs.geodetic_crs, always_xy=True) # Transformer UTM -> lonlat
    

    def lonlat_to_xy(self, lon, lat):
        """
        Computes the corresponding x,y UTM coordinates on the Earth given some lon,lat coordinates.
        This computation is done using the pyproj package.

        Parameters
        ----------
        lon : array ([n_points])
              Values in lon coordinate.
        lat : array ([n_points]) 
              Values in lat coordinate

        Returns
        -------
        x : array ([n_points]) 
            Corresponding values in x coordinate.
        y : array ([n_points])
            Corresponding values in y coordinate.

        """
        (x, y) = self.lonlat_xy_proj.transform(lon, lat)
        return x, y
    
    def xy_to_lonlat(self, x, y):
        """
        Computes the corresponding lon,lat coordinates on the Earth given some x,y UTM coordinates.
        This computation is done using the pyproj package.

        Parameters
        ----------
        x : array ([n_points])
            Values in x coordinate.
        y : array ([n_points])
            Values in y coordinate.

        Returns
        -------
        lon : array ([n_points])
              Corresponding values in lon coordinate.
        lat : array ([n_points])
              Corresponding values in lat coordinate

        """
        (lon, lat) = self.xy_lonlat_proj.transform(x, y)
        return lon, lat

    def assign_lonlat_to_xy(self, dataset):
        """
        Assign new x,y coordinates to dataset using the current lon,lat coordinates

        Parameters
        ----------
        dataset : xarray dataset
                  An xarray dataset containing lon,lat coordinates.

        Returns
        -------
        dataset : xarray dataset 
                  An xarray dataset containing lon,lat coordinates with corresponding x,y coordinates.
        """
        x, y = self.lonlat_to_xy(dataset['lon'].data, dataset['lat'].data)
        dataset = dataset.rename({'lon':'x',
                                  'lat':'y'})
        dataset = dataset.assign_coords({'x':('x', x),
                                         'y':('y', y)})
        return dataset
    
    def assign_xy_to_lonlat(self, dataset):
        """
        Assign new lon,lat coordinates to dataset using the current x,y coordinates

        Parameters
        ----------
        dataset : xarray dataset
                  An xarray dataset containing x,y coordinates.

        Returns
        -------
        dataset : xarray dataset
                  An xarray dataset containing x,y coordinates with corresponding lon,lat coordinates.

        """
        lon, lat = self.xy_to_lonlat(dataset['x'].data, dataset['y'].data)
        dataset = dataset.rename({'x':'lon',
                                  'y':'lat'})
        dataset = dataset.assign_coords({'lon':('lon', lon),
                                         'lat':('lat', lat)})                         
        return dataset
    
