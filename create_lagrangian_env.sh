conda create -y --name lagrangian  python=3.7.0
source activate lagrangian

conda install -y -c conda-forge -c anaconda -c fbriol python-dateutil netcdf4 dask libboost=1.71.0 lagrangian 
conda install xarray matplotlib cartopy
pip install notebook
