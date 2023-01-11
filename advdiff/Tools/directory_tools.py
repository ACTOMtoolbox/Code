# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:36 2022

This module includes tools for creating, cleaning directories and getting files within directories.

@author: Ketil
"""

import glob
import logging
from os     import listdir, remove, mkdir
from shutil import rmtree
from Tools.parser_tools import any_digits

def get_files(data_dir, prefix=None, include=None, ignore=None, filetype='nc', verbose=True):
    """
    Get filesnames from a specific directory, given a prefix to look for, a filetype and what keywords to ignore.

    Parameters
    ----------
    data_dir : str
               Directory path.
    prefix   : str or None, optional
               Prefix for files to look for, by default None. If None, select everything in directory.
    include  : str or None, optional
               Prefix of what filenames must include. If None, select all.
    ignore   : str or None, optional
               Prefix of what filenames to ignore, by default None. If None, do not ignore anything.
    filetype : str, optional
               Filetype to look for, by default 'nc'.

    Returns
    -------
    list or None
                List of filenames found matching the criteria given or None if no files found.
    """
    
    if prefix is None:
        filepath = data_dir + '*.{}'.format(filetype)
    else:
        filepath = data_dir + '*{}*.{}'.format(prefix, filetype)
    if ignore is not None:
        ignore = data_dir + '*{}*.{}'.format(ignore, filetype)
    else:
        ignore = ''
    if verbose:
        logging.info(f'Reading: {filepath}')
    files = sorted(set(glob.glob(filepath)) - set(glob.glob(ignore)))
    if include is not None:
        files = [file for file in files if include in file]
    if not files:
        if verbose:
            logging.info('No files of the form {} found...\n'.format(filepath))
        return None
    else:
        return files


def make_directory(dir_name, verbose=True):
    """
    Make a directory

    Parameters
    ----------
    dir_name : str
               Directory name.
    """
    try:    # Create a folder
        mkdir(dir_name)
    except: # Folder already exists
        if verbose:
            logging.warning(f'The directory {dir_name} already exists.')


def data_cleanup(config, fill_empty=True):
    """
    Delete and clean up outdata, figures and coordinate folders before a simulation is run.

    Parameters
    ----------
    config     : config
                 Config file containing directory paths to the folders of interest.
    fill_empty : bool
                 To fill empty directories with a filler file upon cleanup to ensure directories are stored to github.
    """
    #Delete all files in \Outdata
    file_names=listdir(config['paths']['outdata_path'])
    for file in file_names:
        try:
            rmtree(config['paths']['outdata_path']+file)
        except:
            remove(config['paths']['outdata_path']+file)
    
    #Delete all files in \Figures        
    file_names=listdir(config['paths']['figures_path'])
    for file in file_names:
        try:
            rmtree(config['paths']['figures_path']+file)
        except:
            remove(config['paths']['figures_path']+file)
    
    #Delete all files in \Indata\Coordinates
    in_path = config['paths']['indata_path']+config['paths']['coord_path']
    file_names=listdir(in_path)
    for file in file_names:
        try:
            rmtree(in_path+file)
        except:
            remove(in_path+file)

    #Fill empty directories 
    if fill_empty:
        fill_empty_directories(config=config)


def fill_empty_directories(config):
    """
    Fill folders that become empty upon cleanup. This is to make sure there is some filler text file inside necessary directories.
    This is because empty directories are removed on github.

    Parameters
    ----------
    config : config
             Config file containing directory paths to the folders of interest.
    """
    with open(config['paths']['outdata_path']+'Outdata_goes_here.txt', 'w') as f:
        f.write('Filler file for outdata directory')

    with open(config['paths']['figures_path']+'Figures_go_here.txt', 'w') as f:
        f.write('Filler file for figures directory')

    with open(config['paths']['indata_path']+config['paths']['coord_path']+'UTM_Coordinate_transforms_go_here.txt', 'w') as f:
        f.write('Filler file for coordinate directory')
