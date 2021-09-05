PyUpsc_3D: A Python-based Corner Point Grid Upscaling Library
=========================================================================
Mustapha Zakari(zakari1@univ-lorraine.fr)<br>
*CNRS, Observatoire Terre et environnement de Lorraine, RING TEAM, Nancy*

`PyUpsc_3D` is an open source library for testing upscaling algorithms on cartesian and corner point grids. 
It is based on PyGRDECL for CPG visualization with vtk based visualizers (here pyvista). 
PyGRDECL is a light-weight open source library developped by Bin Wang (binwang.0213@gmail.com), Craft and Hawkins Department of Petroleum Engineering, US.

## Install & Usage
`Anaconda` (https://www.anaconda.com/download/), `vtk`(https://anaconda.org/anaconda/vtk) and `Pyvista` ([https://docs.pyvista.org/getting-started/installation.html](https://docs.pyvista.org/getting-started/installation.html)) is required. The library works on both Windows and Linux.

Following install script can be used for creating a seperate virtual env
```python
#Create Notebook Python environment from python or bash
conda create -n upscale_3D python=3.8
conda activate upscale_3D
conda install jupyter
pip install numpy scipy matplotlib pandas
pip install vtk pyvista itkwidgets
conda install shapely
pip install https://github.com/enthought/mayavi/zipball/master

#For windows users with python3.8 PyUpsc_3D can be directly installed and launched using the "install.bat" and "start.bat" files.
```

# License
PyGRDECL and PyUpsc_3D are released under the terms of the BSD license, and thus free for commercial and research use. Feel free to use these codes into your own project with a PROPER REFERENCE.  

M.Zakari, PyUpsc_3D A Python-based Corner Point Grid Upscaling Library (2021), GitHub repository, https://github.com/mzakari31/PyGRDECL
B. Wang, PyGRDECL A Python-based GRDECL Visualization Library, (2018), GitHub repository, https://github.com/BinWang0213/PyGRDECL

