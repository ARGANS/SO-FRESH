#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
#| James Hickson | Argans UK | jhickson@argans.co.uk |

# Package loader # 
import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import configparser, inspect

# Toolkit loader # 
from data import parser_toolkit as p_tk

import sys,os
parser = argparse.ArgumentParser(description="""
-------------------------------------------------
| Supply the filepath to the configuration file |
-------------------------------------------------
""",
formatter_class=argparse.RawDescriptionHelpFormatter)

"""----------------------------------------------
~~~ Functions ~~~
----------------------------------------------"""
class configuration_reader:
    def __init(self):
        self.check = True

    def config_setup(self):
        """
        Set variables.
        """
        try:
            self.process = config["parameters"]["process"]
            self.data_directory = config["parameters"]["data_directory"]
            self.product = config["parameters"]["product"]
            self.aoi = config["parameters"]["aoi"]
            self.start_date = config["parameters"]["start_date"]
            self.end_date = config["parameters"]["end_date"]
            self.resample = config["parameters"]["resample"]
        except Exception as e:
            raise RuntimeError(f"Configuration exception in '{inspect.stack()[0][3]}': {e}")
    
    def sort_inputs(self, process, data_directory, product, aoi, start_date, end_date, resample):
        """
        Sort inputs and feed to processes.
        """
        # Processes | products | resample
        process = process.split(",")
        product = product.split(",")
        if resample == "False": resample = bool(False)
        elif resample != "False": 
            resample = [float(r) for r in resample.split(",")]
            if len(resample) != 2: raise RuntimeError(f"Please specify two values for an X and Y resample resolution.")
        
        ### ACQUIRE & ADJUST ###
        if any("acquire" in p for p in process) and any("adjust" in p for p in process):
            print("Parsing inputs to 'acquire & adjust'...\n")
            if len(product) > 1: raise RuntimeError(f"'{process[0]} and {process[1]}' only takes one product as an argument, please refine.")
            else: pass
            p_tk.argument_receiver(process, data_directory, product, aoi, start_date, end_date, resample).acquire_adjust_parser()
        ### ACQUIRE ###
        elif any("acquire" in p for p in process): 
            print("Parsing inputs to 'acquire'...\n")
            if len(product) > 1: raise RuntimeError(f"'{process}' only takes one product as an argument, please refine.")
            else: pass
            p_tk.argument_receiver(process, data_directory, product, aoi, start_date, end_date, resample).acquire_parser()
        ### ADJUST ###
        elif any("adjust" in p for p in process): 
            print("Parsing inputs to 'adjust'...\n")
            if len(product) > 1: raise RuntimeError(f"'{process}' only takes one product as an argument, please refine.")
            else: pass
            p_tk.argument_receiver(process, data_directory, product, aoi, start_date, end_date, resample).adjust_parser()
        ### FUSION ###
        elif any("fusion" in p for p in process): 
            print("Parsing inputs to 'fusion'...\n")
            if len(product) < 2: raise RuntimeError(f"'{process}' requires multiple products to fuse, please refine.")
            else: pass
            p_tk.argument_receiver(process, data_directory, product, aoi, start_date, end_date, resample).fusion_parser()
        ### Raise error for impossible combinations. ###
        elif any("acquire" in p for p in process) and any("fusion" in p for p in process):
            raise RuntimeError(f"'{process[0]} and {process[1]}' is not a possible combination.")
        elif any("adjust" in p for p in process) and any("fusion" in p for p in process):
            raise RuntimeError(f"'{process[0]} and {process[1]}' is not a possible combination.")

        '''
        try:
            if process == "acquire":
                product = product[0]
                # (process, data_directory, product, aoi, start_date, end_date)
                # Pass to acquire function toolkit.
            elif process == "adjust":
                product = product[0]
                # (process, data_directory, product, aoi, start_date, end_date, resample)
            elif process == "fusion":
                print()
                # # (process, data_directory, product, aoi, start_date, end_date, resample)
            
        except Exception as e:
            raise RuntimeError(f"The requested process is not possible, please consider refining to any of the following: acquire, adjust or fusion.")
        '''
"""----------------------------------------------
~~~ Main ~~~
----------------------------------------------"""
if __name__ == "__main__":
    try:
        parser.add_argument("config", type=str, help="Full filepath to the configuration file.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
        configfile = args.config

        # Check if configuration file exists.
        if not os.path.isfile(configfile):
            raise RuntimeError("The config file does not exist.")

        # Read configuration file.
        try:
            config = configparser.ConfigParser()
            config.read(configfile)
            print('Configuration file read...\n')
        except Exception as e:
            raiseRuntimeError("Unable to find or read the configuration file. Please review the filpath and inputs.")

        # Set-off processing to appropriate scripts.
        reader = configuration_reader()
        reader.config_setup()
        reader.sort_inputs(reader.process, reader.data_directory, reader.product, reader.aoi, reader.start_date, reader.end_date, reader.resample)

    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)