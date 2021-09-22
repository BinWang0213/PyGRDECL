py -3.8 -m venv env
call env\Scripts\activate
python -m pip install --upgrade pip
python -m pip install jupyter notebook numpy matplotlib 
python -m pip install shapely vtk scipy pyvista itkwidgets
python -m pip install ipyvtklink ipykernel
pause