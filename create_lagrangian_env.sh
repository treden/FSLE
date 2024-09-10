conda create -y --name lagrangian  python=3.7.0
source activate lagrangian

conda install fbriol::lagrangian python-dateutil netcdf4 dask libboost=1.71.0
conda install xarray matplotlib cartopy
pip install notebook
