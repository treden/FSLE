# FSLE Calculation Toolbox

## Introduction

This **FSLE Calculation Toolbox** is a set of overlying functions based on the [Aviso Lagrangian Python Toolbox](https://github.com/CNES/aviso-lagrangian). It help using the functionality of the Aviso-lagrangian toolbox by adding overlaying functions that enable straightforward application with various altimetry datasets and 2D/3D circulation model outputs.

## Features

- Compute 2D and 3D FSLE from various datasets, based on xarray integration.
- Utilizes the Aviso-Lagrangian Python Toolbox for FSLE calculations (https://github.com/CNES/aviso-lagrangian).
- Only regular grid are valid inputs for now (longitude (n) x latitude (m)). The next step should be to manage curvilinear grids.

## Requirements

- Python 3.7
- NumPy
- xarray
- pandas
- [Aviso Lagrangian Python Toolbox](https://github.com/CNES/aviso-lagrangian)

See the **create_lagrangian_env.sh** file to create a working conda environment

## Usage

## Function Descriptions

### `compute_FSLE`

- Compute 2D FSLE for U and V data.
- Specify temporal, spatial, and integration parameters.
- Save results to a NetCDF file.
- Manage temporary files and execute the FSLE calculation.

### `compute_FSLE_3D`

- Compute 3D FSLE for U and V data at multiple depth layers.
- Specify temporal, spatial, and integration parameters.
- Save results in a 3D dataset with depth as a dimension.
- Iterate over depth levels, combine results, and handle temporary files.

## Acknowledgments

- [Aviso Lagrangian Python Toolbox](https://github.com/CNES/aviso-lagrangian)
