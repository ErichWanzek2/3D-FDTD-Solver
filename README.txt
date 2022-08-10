

This software is capable of simulating maxwells equations in a 3D rectangular cavity
with PEC boundaries.
  

Language: Python 

scriptname: FDTD_3D.py

    
Written by: Erich Wanzek
Written 3/26/2021
Intro to Computational Electromangtics
University of Colrado Denver


To run properly, the libraries:

numpy
matplotlib
scipy

must be installed.


This code was developed within and performed within the spyder IDE for python. 
The spyder IDE has the requried libraries already available with its installation.


To input data into the software, there is a txt file that contains all the parameters and input controlls.
This txt file is called FDTD_simulation_parameters.txt
This file must be in the same path folder as the python script FDTD_3D.py so it can be read in by the FDTD_3D.py script.


Here is the layout of the input txt file, with example data:

Lx=1                              #cavity dimenesions in meters
Ly=1
Lz=1  
nx=21                             #number of FDTD nodes along each dimension
ny=21
nz=21
CFLN=0.99                         #CFLN number. Must be less then 1
simtime=4e-7    
source_signature=dgaussian        #source signature type selecton. Defaulted to differentiated gaussian
pulse_delay=5e-9                  #pulse delay of the source signature in seconds
pulse_half_width=1e-9             #pulse half width of the source signature in seconds
modulation_frequency=10e9         #modulation frequency of the source signature in hertz
current_source_number=1           #number of current sources
i1_src=8                          #nodes bounds of source
j1_src=8
k1_src=8
i2_src=8
j2_src=8
k2_src=12
i1_fld=16                         #node bounds of location of field measurement
j1_fld=16
k1_fld=16
i2_fld=16
j2_fld=16
k2_fld=17














 