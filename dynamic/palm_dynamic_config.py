#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------#
#
# Scripts for processing of WRF-CHEM files to PALM dynamic driver.
#
#------------------------------------------------------------------------------#
'''Configuration module.
Configuration options are sourced into this module's globals.
'''
import os.path
from pathlib import Path
import inspect

# Just for PyCharm and similar IDEs to allow autocompletion from config values
if False:
    from palm_dynamic_defaults import *
    from palm_dynamic_init import *

def configure(configname):
    global dir_scripts
    # get path of the palm_dynamic script to source default config and init
    dir_scripts = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    print('\nRunning palm_dynamic from:', dir_scripts)
    # use config defaults
    configdefaultsfile = os.path.join(dir_scripts, "palm_dynamic_defaults.py")
    print('Default case config: ', configdefaultsfile)
    # user config file is located in configurations directory
    configfile = os.path.join(os.path.abspath(Path(dir_scripts).parent), "configurations", configname + '.conf')
    print('User case config:', configfile)
    # initialization of standard parameters done in script
    standardinitfile = os.path.join(dir_scripts, "palm_dynamic_init.py")
    print('Standard initialization: ', standardinitfile)
    # check existence of the supplied config file
    if not os.path.isfile(configfile):
        print("Config file " + configfile + " does not exists!")
        exit(2)
    # read default config values
    exec(open(configdefaultsfile).read(), globals())
    # read user configuration
    exec(open(configfile).read(), globals())
    # perform the standard initialization
    exec(open(standardinitfile).read(), globals())
