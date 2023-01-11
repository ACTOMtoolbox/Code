# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:30:21 2022

This module includes tools for creating artificial sources, probes and velocities.

@author: Ketil
"""

import numpy  as np
import pandas as pd
import xarray as xr
from ast import literal_eval

def create_sources(config):

    coordinate_system      = str(literal_eval(config['setup']['coordinate_system']))
    source_locations       = np.array(literal_eval(config['sources']['locations']))
    location_probabilities = np.array(literal_eval(config['sources']['location_probabilities']))

    Horizontal = source_locations[:,0]
    Vertical   = source_locations[:,1]

    if coordinate_system == 'UTM':
        Hcoord = 'x'
        Vcoord = 'y'
    elif coordinate_system == 'WGS84':
        Hcoord = 'lon'
        Vcoord = 'lat'
    
    data_vars = {'location_probability' : (['source'], location_probabilities)}
    coords    = {Hcoord : (['source'], Horizontal), 
                 Vcoord : (['source'], Vertical)}
    sources = xr.Dataset(data_vars=data_vars, coords=coords)
    attrs = {'Source file' : 'Custom file'}
    sources.attrs = attrs
    return sources

def create_probes(config):

    coordinate_system  = str(literal_eval(config['setup']['coordinate_system']))
    probe_locations    = np.array(literal_eval(config['sources']['locations']))

    Horizontal = probe_locations[:,0]
    Vertical   = probe_locations[:,1]

    if coordinate_system == 'UTM':
        Hcoord = 'x'
        Vcoord = 'y'
    elif coordinate_system == 'WGS84':
        Hcoord = 'lon'
        Vcoord = 'lat'
    
    data_vars = {}
    coords    = {Hcoord : (['source'], Horizontal), 
                 Vcoord : (['source'], Vertical)}
    probes = xr.Dataset(data_vars=data_vars, coords=coords)
    attrs = {'Probe file' : 'Custom file'}
    probes.attrs = attrs

def create_velocity(config):

    coordinate_system  = str(literal_eval(config['setup']['coordinate_system']))
    start_date         = literal_eval(config['velocity']['start_date'])
    end_date           = literal_eval(config['velocity']['end_date'])
    freq               = literal_eval(config['velocity']['time_frequency'])
    Lh                 = literal_eval(config['velocity']['Lh'])
    Lv                 = literal_eval(config['velocity']['Lv'])
    H_mid              = literal_eval(config['velocity']['H_mid'])
    V_mid              = literal_eval(config['velocity']['V_mid'])
    u_vel              = literal_eval(config['velocity']['u_vel'])
    v_vel              = literal_eval(config['velocity']['v_vel'])

    if coordinate_system == 'UTM':
        Hcoord = 'x'
        Vcoord = 'y'
    elif coordinate_system == 'WGS84':
        Hcoord = 'lon'
        Vcoord = 'lat'

    time = pd.date_range(start=np.datetime64(start_date), end=np.datetime64(end_date), freq=freq)

    N_nodes = 100
    H = Lh*np.linspace(0,1,int(np.sqrt(N_nodes))) - Lh/2 + H_mid
    V = Lv*np.linspace(0,1,int(np.sqrt(N_nodes))) - Lv/2 + V_mid

    H, V = np.meshgrid(H, V)
    Horizontal = H.flatten()
    Vertical   = V.flatten()

    u = u_vel * np.ones(shape=(time.shape[0],  H.shape[0]))
    v = v_vel * np.ones(shape=(time.shape[0],  V.shape[0]))

    data_vars = {'u' : (['time','node'], u),
                 'v' : (['time','node'], v)}

    velocity = xr.Dataset(data_vars = data_vars, 
                          coords    = {Hcoord    : (['node'], Horizontal), 
                                       Vcoord    : (['node'], Vertical), 
                                       "time"    : (['time'], time),})

    attrs = {'Velocity file' : 'Custom file'}
    velocity.attrs = attrs
    return velocity



