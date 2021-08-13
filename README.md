# geodyn_cmbtopo_models

Python scripts for processing and vizualizing core-mantle boundary topography models.

# Installation
1. Clone this directory
2. Using [conda](https://anaconda.org/), create the ```geotopo``` environement with the required dependencies
```
conda env create -f environment.yml
conda activate geotopo
```

# Content
- ``rotate_models.py``: rotate spherical harmonics CMB topography models to align degree 2 with S20RTS

# Usage
From the ``topo`` subdirectory in the cloned project folder:
```
conda activate geotopo
python rotate_models.py
```
