import os 
import subprocess
import xarray as xr
import pandas as pd
import numpy as np

def time_to_julian_date(t):
    """
    Convert a pandas Timestamp to Julian date relative to '1950-01-01'.
    """
    return (t).to_julian_date() - pd.Timestamp('1950-01-01').to_julian_date()

def set_lagrangian_attributes(ds, t):
    """
    Set needed variables and attributes in the dataset.
    """
    
    ds = ds.assign_coords(LatLon = np.array([0., 1.]))
    ds = ds.assign(LatLonMin = ('LatLon', np.hstack([ds.latitude.min(), ds.longitude.min()])))
    ds = ds.assign(LatLonStep = ('LatLon', np.hstack([ds.latitude.diff(dim = 'lat')[0], ds.longitude.diff(dim = 'lon')[0]])))

    attribute_dicts = {
            'latitude': {'standard_name': 'latitude', 'long_name': 'Latitude', 'units': 'degrees_north', 'axis': 'Y', 'bounds': 'latitude_bnds'},
            'longitude': {'standard_name': 'longitude', 'long_name': 'Longitude', 'units': 'degrees_east', 'axis': 'X', 'bounds': 'longitude_bnds'},
            'LatLon': {'long_name': 'No sense but necessary for some automatic tools', 'units': 'count'},
            'LatLonMin': {'long_name': 'Latitude/Longitude of south/ouest corner', 'units': 'degree'},
            'LatLonStep': {'long_name': 'latitude/longitude steps', 'units': 'degree'},
            'u': {'long_name': 'U', 'units': 'm/s'},
            'v': {'long_name': 'V', 'units': 'm/s'},
        }
    
    for var, attributes in attribute_dicts.items():
        ds[var].attrs = attributes

    # Add date-specific attributes
    JD = time_to_julian_date(t)
    ds['u'].attrs['Date_CNES_JD'] = JD
    ds['u'].attrs['date'] = str(t)
    ds['v'].attrs['Date_CNES_JD'] = JD
    ds['v'].attrs['date'] = str(t)

    return ds

def create_input_dataset(ds, var, domain=None):
    """
    Create a dataset for input variables with optional domain selection.
    """
    data_vars = {
        'longitude': ds[var['longitude']].values,
        'latitude': ds[var['latitude']].values,
        'u': (["time", "latitude", "longitude"], ds[var['u']].values),
        'v': (["time", "latitude", "longitude"], ds[var['v']].values),
        'time': ds[var['time']]
    }

    ds = xr.Dataset(data_vars)

    if not isinstance(domain, type(None)):
        ds = ds.sel(longitude=slice(domain[0], domain[1]))
        ds = ds.sel(latitude=slice(domain[2], domain[3]))

    return ds

def write_list_ini(uv_files, temp_dir):
    """
    Write a list of U and V file paths to 'list.ini' in the specified temporary directory.

    Parameters:
    - uv_files (list): A list of file paths containing U and V data files.
    - temp_dir (str): The directory where 'list.ini' will be created.

    Returns:
    - None
    """
    with open(f'{temp_dir}/list.ini', 'w+') as f:
        for file in uv_files:
            f.write(f'U = {file}\n')
        for file in uv_files:
            f.write(f'V = {file}\n')
        f.write('U_NAME = u\nV_NAME = v\nFILL_VALUE = 0')
        
        
def create_temp_dir(temp_dir):
    """
    Create a temporary directory if it doesn't exist.

    Parameters:
    - temp_dir (str): The path to the temporary directory.

    Returns:
    - None
    """
    try:
        os.mkdir(temp_dir)
    except OSError as error:
    	pass
#        print(temp_dir, 'already exists')
        
def write_input_files(ds, temp_dir):
    """
    Write individual U and V NetCDF files and 'list.ini' in the specified temporary directory.

    Parameters:
    - ds (xarray.Dataset): The dataset containing U and V data for multiple time steps.
    - temp_dir (str): The directory where the files will be written.

    Returns:
    - None
    """
    uv_files = []
    for i in range(len(ds.time)):
        ds_out = ds.isel(time=i)
        t = pd.Timestamp(ds_out.time.values)
        label = ''.join(str(t.date()).split('-'))
        ds_out = set_lagrangian_attributes(ds_out, t)
        ds_out.to_netcdf(f'{temp_dir}/uv_{label}.nc')
        uv_files.append(f'{temp_dir}/uv_{label}.nc')

    write_list_ini(uv_files, temp_dir)

def compute_FSLE(ds, variables, t0, t1, domain, resolution=0.05, stencil = 'triplet', initial_separation = 0.02, final_separation=0.6, output_file='FSLE.nc', temp_dir='_tmp/', mode= 'fsle'):
    """
    Compute Finite-Size Lyapunov Exponents (FSLE) for U and V data.

    Parameters:
    - ds (xarray.Dataset): The dataset containing U and V data over time.
    - variables (dict): Dictionary specifying variable names for 'time', 'longitude', 'latitude', 'u', and 'v'.
    - t0 (str or pandas.Timestamp): The start time for the analysis.
    - t1 (str or pandas.Timestamp): The end time for the analysis.
    - domain (list or None): Spatial domain in the format [x_min, x_max, y_min, y_max]. None for the full domain.
    - resolution (float): Spatial resolution in degrees.
    - final_separation (float): Separation threshold for FSLE calculation.
    - output_file (str): Output NetCDF file to store the FSLE results.
    - temp_dir (str): Temporary directory for intermediate files.

    Returns:
    - None
    """
    # dt = integration_time
    # t0, t1 = time - pd.Timedelta(days=int(dt / 2 + 1)), time + pd.Timedelta(days=int(dt / 2 + 1))
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    
    ds = create_input_dataset(ds, variables, domain=domain)
    ds = ds.sel(time=slice(t0, t1))

    create_temp_dir(temp_dir)
    write_input_files(ds, temp_dir)

    advection_time = int(((ds.time[-1] - ds.time[0])).values.astype(float) / 1e9 / 86400) - 1

    xmin = ds.longitude.values.min()
    xmax = ds.longitude.values.max()
    ymin = ds.latitude.values.min()
    ymax = ds.latitude.values.max()

    t_init = ''.join(str(pd.Timestamp(ds.time.values[-1]).date()).split('-'))

    FSLE_cmd = f'map_of_fle {temp_dir}list.ini {output_file} --mode {mode}\
        {t_init} --advection_time {advection_time} --resolution {resolution} \
        --stencil {stencil} --x_min {xmin} --x_max {xmax} --y_min {ymin} --y_max {ymax} \
        --initial_separation {initial_separation} --final_separation {final_separation} --time_direction backward'
        
    print(FSLE_cmd)

    subprocess.call(FSLE_cmd, shell=True)
    
def compute_FSLE_3D(ds, variables, t0 , t1, domain, depths=None, resolution=0.05, stencil = 'triplet', initial_separation = 0.02, final_separation=0.6, temp_dir='_tmp/'):
    """
    Compute Finite-Size Lyapunov Exponents (FSLE) for U and V data at multiple depth layers.

    Parameters:
    - ds (xarray.Dataset): The dataset containing U and V data at various depth levels.
    - variables (dict): Dictionary specifying variable names for 'time', 'longitude', 'latitude', 'u', 'v', and 'depth'.
    - t0 (str or pandas.Timestamp): The start time for the analysis.
    - t1 (str or pandas.Timestamp): The end time for the analysis.
    - domain (list or None): Spatial domain in the format [x_min, x_max, y_min, y_max]. None for the full domain.
    - depths (list or None): List of depth levels to compute FSLE for. If None, it uses all available depths.
    - resolution (float): Spatial resolution in degrees.
    - final_separation (float): Separation threshold for FSLE calculation.
    - temp_dir (str): Temporary directory for intermediate files.

    Returns:
    - fsle_3d (xarray.Dataset): 3D dataset with depth as a dimension and 'lambda1' variable containing FSLE results.
    """
    fsle_list = []

    if depths is None:
        depths = ds[variables['depth']].values

    for depth in depths:
        output_file = f'FSLE_depth_{depth}.nc'
        compute_FSLE(ds.interp(depth=depth, method='nearest'), variables, t0, t1, domain, resolution, stencil, initial_separation, final_separation, output_file, temp_dir)

        fsle_layer = xr.open_dataset(output_file)
        fsle_list.append(fsle_layer)
        os.remove(output_file)

    fsle_3d = xr.concat(fsle_list, dim=xr.DataArray(depths, name='depth', dims='depth'))
    return fsle_3d
