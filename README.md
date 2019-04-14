PyGRDECL: A Python-based GRDECL Visualization Library
==============================================================================================
Bin Wang (binwang.0213@gmail.com), Craft and Hawkins Department of Petroleum Engineering, Louisiana State University, Baton Rouge, US

<p align="center">
  <img src = "https://github.com/BinWang0213/PyGRDECL/blob/master/img/GridPreview.png" height="300">
</p>

`PyGRDECL` is a light-weight open source library for converting a Eclipse/Petrel grid with properties to a vtu-file.
(to be opened in ParaView for example). 

## Install & Usage

`Anaconda 5.3` (https://www.anaconda.com/download/), `vtk`(https://anaconda.org/anaconda/vtk) and `Shapely 1.5` ([https://anaconda.org/scitools/shapely](https://anaconda.org/conda-forge/shapely)) is required. The library works on both Windows and Linux.

After downloading and unzipping the current <a href="https://github.com/BinWang0213/PyGRDECL/archive/master.zip">repository</a>, navigate to the library directory and run a simple example contained in `Example_GettingStart.ipynb`:

Following install script can be used for creating a seperate virtual env
```python
#Windows
conda create -n pyGRDECL
activate pyGRDECL
conda install pandas numpy matplotlib jupyter notebook scipy 
conda install shapely vtk

#Linux
conda create -n pyGRDECL
source activate pyGRDECL
conda install pandas numpy matplotlib jupyter notebook scipy 
conda install shapely vtk
```

## Getting Start

<p align="center">
  <img src = "https://github.com/BinWang0213/PyGRDECL/blob/master/img/DomeModel.png" height="300">
</p>

```python
from GRDECL2VTK import * 

#Read GRDECL File
Model=GeologyModel(filename='./ExampleData/dome.grdecl')

#Convert ECLIPSE grdecl format into VTK
Model.GRDECL2VTK()

#Decompose the model into sub-volumes in terms of faults automatically (this function requires shapely library)
Model.decomposeModel()

#Output to VTK format
Model.Write2VTU()

#Load a custom new keyword from file
TempData=Model.LoadCellData(varname="TEMP",filename='./ExampleData/dome_Temperature.txt')

#Update model and output to VTK format
Model.Update()
Model.Write2VTU()

```


# License

This code is released under the terms of the BSD license, and thus free for commercial and research use. Feel free to use the code into your own project with a PROPER REFERENCE.  

B. Wang, PyGRDECL A Python-based GRDECL Visualization Library, (2018), GitHub repository, https://github.com/BinWang0213/PyGRDECL

