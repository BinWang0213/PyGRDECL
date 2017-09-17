PyGRDECL: A Python-based GRDECL Visulization Library
==============================================================================================
Bin Wang (binwang.0213@gmail.com), Yin Feng
Department of Petroleum Engineering, Univeristy of Louisiana at Lafayette, Lafayette, US, 70506

`PyGRDECL` is a open source library for converting a Eclipse grid with properties to a vtu-file.
(to be opened in ParaView for example)

# Read a simple grid file

<p align="center">
  <img src = "https://github.com/BinWang0213/PyGRDECL/blob/master/img/PermX.png" height="300">
</p>

After downloading and unzipping the current <a href="https://github.com/BinWang0213/PyGRDECL/archive/master.zip">repository</a>, navigate to the library directory and run a simple example contained in `Example.ipynb`:

<p align="center">
  <img src = "https://github.com/BinWang0213/PyGRDECL/blob/master/img/Fault.png" height="300">
</p>

```python
from GRDECL2VTK import * 

#Read GRDECL File
Grid1=GRDECL_Viewer(filename='./Example/HW1.GRDECL',nx=15,ny=8,nz=1)

#Compute transmissibility in x,y,z direction
Grid1.calc_Trans()

#Output to VTK format
Grid1.write_VTU(filename='./Example/HW1_Perm_Poro',mode=0)
Grid1.write_VTU(filename='./Example/HW1_Trans',mode=1)
Grid1.write_VTU(filename='./Example/HW1_Faults',mode=2)
```


# License

This code is released under the terms of the BSD license, and thus free for commercial and research use. Feel free to use the code into your own project with a PROPER REFERENCE.  
