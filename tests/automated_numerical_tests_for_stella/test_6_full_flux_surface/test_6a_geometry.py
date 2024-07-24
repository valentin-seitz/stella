
################################################################################
#                     Check the evolution of the potential                     #
################################################################################
# These test are designed to test the evolution of the potential, starting
# with a linear and nonlinear simulation using Miller geometry.
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                  Check whether the diagnostics are working                   #
#-------------------------------------------------------------------------------
def test_miller_geometry_for_FFS(tmp_path):    

    # Input file name 
    input_filename = 'miller_linear.in'  
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    
    # Check whether the potential data matches in the netcdf file
    keys = ['t', 'phi2']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False)
    print(f'\n  -->  The potential diagnostics are working.')
    return
    
    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.final_fields' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.final_fields'  
     
    # Check whether the txt files match  
    compare_local_txt_with_expected_txt(local_file, expected_file, name='Final fields', error=False)

    # If we made it here the test was run correctly 
    print(f'  -->  The final fields diagnostics are working.')
    return 
    
