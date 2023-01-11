"""
Created on Thu Jun 23 11:33:46 2022

@author: Ketil

Module for runnnign main.py in an ensemble.
"""

import sys
import itertools
import configparser
import logging
from main                     import main
from Tools.directory_tools    import data_cleanup
from Tools.parser_tools       import strtolist
from datetime                 import datetime


if __name__ == '__main__':

    ### READ MASTER CONFIG FILE ###
    masterconfig=configparser.ConfigParser(allow_no_value=True)
    masterconfig.optionxform = str  # preserve case for letters
    masterconfig.read('Indata/setup.ini')

    ### GET A COPY CONFIG FILE TO BE USED FOR EACH SCENARIO ###
    config=configparser.ConfigParser(allow_no_value=True)
    config.optionxform = str  # preserve case for letters
    config.read('Indata/setup.ini')
    
    ### DATA CLEANUP OF PREVIOUS RUNS ###
    data_cleanup(config=masterconfig) # Make sure this comes before anything is stored

    ### CREATE COPY OF CURRENT MASTER CONFIG FILE IN OUTPUT FOR REPRODUCIBILITY ###
    with open(masterconfig['paths']['outdata_path']+f'master_setup.ini', 'w') as configfile:
        masterconfig.write(configfile)

    ### SET UP LOGGER ###
    levels = {'DEBUG': logging.DEBUG, 'INFO': logging.INFO, 'WARNING': logging.WARNING, 'CRITICAL': logging.CRITICAL}
    logging.basicConfig(format   = '[%(levelname)s] %(message)s', 
                        level    = levels[masterconfig['DEFAULT']['verbose'].upper()], 
                        handlers = [logging.FileHandler(masterconfig['paths']['outdata_path']+'debug.log'), logging.StreamHandler(sys.stdout)])
    logging.info(f'AdvDiff Log: {datetime.now().strftime("%d.%m.%Y-%H.%M.%S")}\n')

    ### SCENARIO VARIABLES ###
    velocity_files     = strtolist(masterconfig['velocity']['velocity_file'], dtype=str)
    time_starts        = strtolist(masterconfig['setup']['time_start'],       dtype=str) 
    time_deltas        = strtolist(masterconfig['setup']['time_delta'],       dtype=str)        
    time_seeds         = strtolist(masterconfig['setup']['time_seed'],        dtype=int)
    BC_types           = strtolist(masterconfig['setup']['BC_type'],          dtype=str)
    Ds                 = strtolist(masterconfig['setup']['D'],                dtype=float)
    diffusion_type     = strtolist(masterconfig['setup']['diffusion_type'],   dtype=str)
    convection_type    = strtolist(masterconfig['setup']['convection_type'],  dtype=str)
    grid_type          = strtolist(masterconfig['grid']['type'],              dtype=str)
    grid_type_out      = strtolist(masterconfig['grid']['type_out'],          dtype=str)
    scenario_variables = [velocity_files, time_starts, time_deltas, time_seeds, BC_types, Ds, diffusion_type, convection_type, grid_type, grid_type_out]
    scenarios          = list(itertools.product(*scenario_variables))

    
    for ii, scenario in enumerate(scenarios):
        logging.info(f'---- ### STARTING SCENATIO: {ii} ### ----\n')

        suffix = f'Scenario_{ii}/' # Name tag of sub-folder to store outputs in. If an empty string, then store outputs in main outdata folder

        velocity_file   = scenario[0]
        time_start      = scenario[1]
        time_delta      = scenario[2]
        time_seed       = scenario[3]
        BC_type         = scenario[4]
        Diff            = scenario[5]
        diffusion_type  = scenario[6]
        convection_type = scenario[7]
        grid_type       = scenario[8]
        grid_type_out   = scenario[9]

        # Set scenario variables in config copy.
        config.set('velocity', 'velocity_file',   str(velocity_file))
        config.set('setup',    'time_start',      str(time_start))
        config.set('setup',    'time_delta',      str(time_delta))
        config.set('setup',    'time_seed',       str(time_seed))
        config.set('setup',    'BC_type',         str(BC_type))
        config.set('setup',    'D',               str(Diff))
        config.set('setup',    'diffusion_type',  str(diffusion_type))
        config.set('setup',    'convection_type', str(convection_type))
        config.set('grid',     'type',            str(grid_type))
        config.set('grid',     'type_out',        str(grid_type_out))
        config.set('paths',    'outdata_path',    masterconfig['paths']['outdata_path']+suffix)
        config.set('paths',    'figures_path',    masterconfig['paths']['figures_path']+suffix)

        ### RUN ADVDIFF MODULE ###
        main(config=config)

        logging.info(f'---- ### SCENARIO: {ii} COMPLETED ### ---\n')
    
    logging.info(f'---- ### ALL SCENARIOS COMPLETED ### ----')

    ### CLOSE LOGGER ###
    logging.shutdown()