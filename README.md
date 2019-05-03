# codes-dynamic-rupture-and-earthquakes
Codes developped by Simon Guerin-Marthe during the PhD of Geophysics at Durham University

list of codes and brief description:

-Main-plotter.py and read_function.py : 

Main-plotter.py calls read_function.py, 
which together enable to read all the datasets from hdf5 files produces by the
high-speed acquisition system. Used mainly to process strain-gage signals.

-Iris_EQ_animation:

Python script (jupyter notebook) which finds online earthquake datasets from Iris website,
and generates animations of earthquakes sequences for chosen timespans and regions.

-dynamic_fd_intro.m

Matlab script which simulates dynamic inplane ruptures

-dynamic_rupture_snapshot.ipynb

Python script (jupyter notebook) for the simulation of dynamic inplane ruptures.

-DG_version1_static_elasticity_FEM.ipynb

Python jupyter notebook used for static finite element simulations 
with simple quadrilateral elements, and possibility to introduce a crack in the domain.

-FEM2D_models_of_biaxial_loading.ipynb

Python jupyter notebook used for static finite element simulations , and more specifically
for simulation of isochromatics in polycarbonate plates with different loading configurations.

-DGv2_spyder.py

Python script under developement, to simulate dynamic rupture by prescribing tractions at split nodes
along an interface.

-inversion-isochromatics-fringes.ipynb

Python jupyter notebook used to calculate strains in 2D from PIV analysis results 
(text files such as pivImJ.txt) from the PIV plugin of imageJ free software.


