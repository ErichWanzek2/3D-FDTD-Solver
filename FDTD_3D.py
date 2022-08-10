
"""
Created on Wed Feb 17 20:31:45 2021

Maxwell's equations simulator for rectangular PEC cavity using 3D FDTD Yee 
Method

    This software is capable of simulating Maxwells euations in a three 
    dimensional rectanguler PEC cavity using the FDTD Yee Method. The software
    can also calcualte numericaly the resonsant modes of the cavity and comapre
    them to the anayltical resonsant modes of the theoritecal PEC rectangular
    cavity.
     
Written by: Erich Wanzek
Written 3/29/2021
Intro to Computational Electromangtics
University of Colrado Denver
"""
##############################################################################
##############################################################################
##############################################################################
'''Import Required Libraries'''
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
##############################################################################
plt.rcParams['figure.dpi']=500
plt.style.use('default')
##############################################################################
##############################################################################
def TimeAdvance():
    '''This function time-advances the main electric and magnetic field nodes
    and also injects the source electreic field, and save the field data at 
    the specified nodes.'''
    
    for t in range(0,nt):
        'Update electric fields throughout entire simulation volume, EXCCLUDING boundaries'
        E_Field_Update();
        
        'Electric field source injection.'
        E_Field_Source_Injection(t);
        
        'Update magenetic fields throughout entire simulation volume'
        H_Field_Update();
        
        'Save output data'
        save_output(t); 
        
    return None
##############################################################################
def E_Field_Update():
    """This function updates the main E field vector components at the simulation
    nodes by using the FDTD yee method update equations for Ampere-Maxwell
    equaiton"""
    
    '''Main update loops for the electric field vector components'''
    for k in range(1,nz-1):
        for j in range(1,ny-1):
            for i in range(0,nx-1):
                Ex[i][j][k] = Ex[i][j][k] + cExy*(Hz[i][j][k] - Hz[i][j-1][k]) - cExz*(Hy[i][j][k] - Hy[i][j][k-1]);
    
    for k in range(1,nz-1):
        for j in range(0,ny-1):
            for i in range(1,nx-1):
                Ey[i][j][k] = Ey[i][j][k] + cEyz*(Hx[i][j][k] - Hx[i][j][k-1]) - cEyx*(Hz[i][j][k] - Hz[i-1][j][k]);
    
    for k in range(0,nz-1):
        for j in range(1,ny-1):
            for i in range(1,nx-1):
                Ez[i][j][k] = Ez[i][j][k] + cEzx*(Hy[i][j][k] - Hy[i-1][j][k]) - cEzy*(Hx[i][j][k] - Hx[i][j-1][k]);
                ##print(Ez[i][j][k])
    return None                
##############################################################################
def H_Field_Update():
    """This function updates the main H field vector components at the simulation
    nodes by using the FDTD yee method update equations for Faraday-Maxwell
    equaiton"""
    
    '''Main update loops for the magnetic field vector components'''
    for k in range(0,nz-1):
        for j in range(0,ny-1):
            for i in range(0,nx): 
                Hx[i][j][k] = Hx[i][j][k] - cHxy*(Ez[i][j+1][k] - Ez[i][j][k]) + cHxz*(Ey[i][j][k+1] - Ey[i][j][k]);
    
    for k in range(0,nz-1):
        for j in range(0,ny):
            for i in range(0,nx-1):
                Hy[i][j][k] = Hy[i][j][k] - cHyz*(Ex[i][j][k+1] - Ex[i][j][k]) + cHyx*(Ez[i+1][j][k] - Ez[i][j][k]);
    
    for k in range(0,nz):
        for j in range(0,ny-1):
            for i in range(0,nx-1):
                Hz[i][j][k] = Hz[i][j][k] - cHzx*(Ey[i+1][j][k] - Ey[i][j][k]) + cHzy*(Ex[i][j+1][k] - Ex[i][j][k]);
    
    return None                
##############################################################################
def E_Field_Source_Injection(n):
    """This function injects z vector-directed current source within the 
    bounds specified in the input parameter file. Points are labeled as 
    (i1_src, j1_src, k1_src) to (i2_src, j2_src, k2_src). In the specific
    case of this function, the injected source is a differentiated gaussian
    pulse.
    """
    t=(n-0.5)*dt
    
    Jz=(-2/tw)*(t-t0)*np.exp((-(t-t0)**2)/(tw**2));
    
    for k in range(k1_src-1,k2_src):
        for j in range(j1_src-1,j2_src):
            for i in range(i1_src-1,i2_src):
               Ez[i][j][k] = Ez[i][j][k] + (dt)*Jz;
               #print(Ez[i][j][k])
    return None 
##############################################################################
def save_output(n):
    """This function samples and save the summed z vector-directed total 
    electric field within the bounds specified in the input parameter file. 
    The sample volume span is bounded by the points labeled as 
    (i1_fld, j1_fld, k1_fld) to (i2_fld, j2_fld, k2_fld)."""
    
    for k in range(k1_fld-1,k2_fld):
        for j in range(j1_fld-1,j2_fld):
            for i in range(i1_fld-1,i2_fld):
               Ez_out[0][n] = Ez_out[0][n] +  Ez[i][j][k]
               #print(Ez_out)
    return None 
    
##############################################################################
def readin_parameter_data_file():
    """The function readin_parameter_data_file() is a function that reads in 
    the inital simualtion variables from the text file called
    FDTD_simulation_parameters.txt. The software user must edit the values 
    of the variables in this txt file and save the file before running the 
    python simulation script. The function readin_parameter_data_file() is 
    called once at the beginning of the main() function"""
    
    variables={}
    with open('FDTD_simulation_parameters.txt') as f:
         for line in f:
             name, value = line.strip().split('=')
           
             if name == 'source_signature':
                 variables[name]= value
             elif name == 'make_video':
                 variables[name]= value
             elif name == 'graph_source_voltage':
                 variables[name]= value 
             elif name == 'graph_load_voltage':
                 variables[name]= value
             elif name == 'graph_line_voltage':
                 variables[name]= value     
             else:    
                 variables[name]= float(value)
             
    global Lx; Lx = variables['Lx']       #Global Cavity Dimesnions a, b, c
    global Ly; Ly = variables['Ly']
    global Lz; Lz = variables['Lz']
    
    global nx; nx = int(variables['nx'])       #Discrete sampling Density:Number of nodes in each dimension span
    global ny; ny = int(variables['ny'])
    global nz; nz = int(variables['nz'])
    
    global CFLN; CFLN = variables['CFLN']           #courant number
    global simtime; simtime = variables['simtime']  #simulation time
    

    global source_sig; source_sig= variables['source_signature']
    
    global t0; t0 = variables['pulse_delay']
    global tw; tw = variables['pulse_half_width']
    global modulation_frequency; modulation_frequency= variables['modulation_frequency']
    
    
    global current_source_number; current_source_number=variables['current_source_number']

    global i1_src; i1_src = int(variables['i1_src'])
    global j1_src; j1_src = int(variables['j1_src'])
    global k1_src; k1_src = int(variables['k1_src'])
    
    global i2_src; i2_src = int(variables['i2_src'])
    global j2_src; j2_src = int(variables['j2_src'])
    global k2_src; k2_src = int(variables['k2_src'])
    
    
    global i1_fld; i1_fld = int(variables['i1_fld'])
    global j1_fld; j1_fld = int(variables['j1_fld'])
    global k1_fld; k1_fld = int(variables['k1_fld'])
    
    global i2_fld; i2_fld = int(variables['i2_fld'])
    global j2_fld; j2_fld = int(variables['j2_fld'])
    global k2_fld; k2_fld = int(variables['k2_fld'])

    return None
##############################################################################
def initialize_parameters():
    """The function initialize_parameters() precalculates the physcial 
    constants and FDTD node/grid parameters from the parameters inputted in 
    the input data file.This funciton is performed once when 
    initialize_parameters() is called in setup_FDTD_sim()"""
    
    global c; c=2.99782458*10**8;  #speed of light m/s
    global etao; etao=376.7303134617706554679; #impedance of free space Ohms
    
    global muo; muo= etao/c;
    global epso; epso=1/(etao*c);
    
    
    'compute grid discretization parameters'
    global dx; dx = Lx/(nx-1);
    global dy; dy = Ly/(ny-1);
    global dz; dz = Lz/(nz-1);
    
    global dt; dt= CFLN /(c*(np.sqrt((1/(dx**2))+(1/(dy**2))+(1/(dz**2)))))
    global nt; nt= m.floor(simtime/dt);
    print(nt)
    
    return None

##############################################################################
def initialize_arrays():
    """The function initialize_arrays() initializes and sets up the data 
    arrays for the main E field and main H Field vector component 
    values. The function also initializes the time and output arrays.
    This funciton is performed once when initialize_arrays() is called 
    in setup_FDTD_sim()"""
    
    'Initialize the E_Field and H_Field vector component arrays'
    global Ex; Ex=np.zeros((nx-1,ny,nz));
    global Ey; Ey=np.zeros((nx,ny-1,nz));
    global Ez; Ez=np.zeros((nx,ny,nz-1));
 
    global Hx; Hx=np.zeros((nx,ny-1,nz-1));
    global Hy; Hy=np.zeros((nx-1,ny,nz-1));
    global Hz; Hz=np.zeros((nx-1,ny-1,nz));
    print(Hx.shape)
    "Initialize Ez out for ouptu data storage for postprocessing"
    global Ez_out; Ez_out=np.zeros((1,nt));

    '''Initialize  time array'''
    global time; time=np.zeros((1,nt))
   
    '''Populate distance and time arrays with space/time values'''
    for t in range(0,nt):    
        time[0][t]=(dt*(t))
    
    return None
##############################################################################
def calculate_prefactors():
    """The function calculate_prefactors() calcualtes the prefactor 
    coeffiecents for the FDTD update equations for the field vecotr compents.
    This funciton is performed once when calculate_prefactors() is called 
    in setup_FDTD_sim()"""
    
    '''calcualte the prefactor coefficeint used in the FDTD update equationns'''
    global cExy; cExy=(dt/(epso*dy));
    global cExz; cExz=(dt/(epso*dz));
    global cEyz; cEyz=(dt/(epso*dz));
    global cEyx; cEyx=(dt/(epso*dx));
    global cEzx; cEzx=(dt/(epso*dx));
    global cEzy; cEzy=(dt/(epso*dy));
    
    global cHxy; cHxy=(dt/(muo*dy));
    global cHxz; cHxz=(dt/(muo*dz));
    global cHyz; cHyz=(dt/(muo*dz));
    global cHyx; cHyx=(dt/(muo*dx));
    global cHzx; cHzx=(dt/(muo*dx));
    global cHzy; cHzy=(dt/(muo*dy));
     
    return None
##############################################################################
def setup_FDTD_sim():
    """setup_FDTD_sim function is the main setup function that runs once in 
    the main funciton, when main() is called. This function runs the setup 
    subroutine functions to initialize the data arrays, calcualte simulation 
    parameters, and calcualte FDTD update equation prefactors"""
    
    initialize_parameters();
    initialize_arrays();
    calculate_prefactors();

    return None
##############################################################################
##############################################################################
"""Auxilary functions for post processing"""
##############################################################################
def Fourier_Specturm():
    """The function Plot_Electric_Field() plots the summed total z-component
    transient electric field in the node volume boundend by the nodes labeld
    as:(i1_fld, j1_fld, k1_fld) to (i2_fld, j2_fld, k2_fld) which are 
    specified in the input parameter txt file."""
    
    Nfft=2**12 #2**18 ##nt*8;
    print(Nfft)
    df=1/(Nfft*dt);
    
    FT_Ez=np.absolute(np.fft.fft(Ez_out[0],Nfft))
    print(len(FT_Ez))
    fmax=600e6;
    freq_samps=int(np.floor(fmax/df))

    '''Initialize  frequency array'''
    global freq; freq=np.zeros((1,freq_samps))
    '''Populate distance and time arrays with space/time values'''
    for f in range(0,freq_samps):    
        freq[0][f]=(df*(f))
    
        
    order=4
    freq_modes=np.zeros((order,order,order)); 
    
    for mm in range(0,order-1):
        for n in range(0,order-1):
            for l in range(0,order-1):
                freq_modes[mm][n][l]= (c/(2*np.pi)) * (np.sqrt((((mm*np.pi)/Lx)**2)+(((n*np.pi)/Ly)**2)+(((l*np.pi)/Lz)**2)))
                #print(freq_modes[mm][n][l])
         
        
    title=('Fourier Spectrum of Transient Ez Electric Field, Node (16,16,16)-(16,16,17)')
    y_max=1.1*np.amax(FT_Ez[0:freq_samps])   
    fig, ax =plt.subplots()
    ax.plot((freq[0])*1e-6,FT_Ez[0:freq_samps],color='k',lw=0.5)
    ax.grid(linestyle='--')
    ax.set_ylim((0,y_max))
    #ax.set_ylim(-250,-150) # 10*np.log(y_max))
    ax.set_xlim(0,fmax*1e-6)
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Electric Field (V/m)')
    
    print(len(freq[0]))
    #fmaxer=int(np.floor(fmax/df))
    samplecut=int(np.floor((4.5/8)*freq_samps))
    print(samplecut)
    peaks,_ = find_peaks(FT_Ez[0:samplecut],height=0.8e-8,distance=10)
    #peaks,_ = find_peaks(FT_Ez[0:samplecut],height=0.25e-8,distance=100)
    #peaks,_ = find_peaks(FT_Ez[0:samplecut],height=0.02e-8,distance=1000)
    print(peaks)
    
    colors=['g','c','b','r','b','m','k']
    i=0
    for j in peaks:
        frequency=freq[0][j]*1e-6
        ax.plot(frequency, FT_Ez[j], "x", label='%.2f' %frequency +' MHz',color=colors[i])
        i=i+1

    '''
    
    i=0
    for mm in range(0,order-1):
        for n in range(0,order-1):
            for l in range(0,order-1):
                if mm ==0:
                    if n==0: 
                        n=n+1
                    if l==0: 
                        l=l+1
                if n ==0:
                    if mm==0:
                        mm=mm+1
                    if l==0: 
                        l=l+1
                if l ==0:
                    if mm==0:
                        mm=mm+1
                    if n==0: 
                        n=n+1
                print(mm,n,l)
                ax.axvline(x=freq_modes[mm][n][l]*1e-6,lw=0.7,color=colors[i], label=(str(mm)+','+str(n)+','+str(l)))
                i=i+1
    '''  
            
    ax.axvline(x=freq_modes[1][1][0]*1e-6,lw=0.6,color='g',label='f 110 = ' + '%.2f' %(freq_modes[1][1][0]*1e-6) +' MHz' )
    ax.axvline(x=freq_modes[1][1][1]*1e-6,lw=0.6,color='c',label='f 111 = ' + '%.2f' %(freq_modes[1][1][1]*1e-6) +' MHz' )
    ax.axvline(x=freq_modes[2][1][0]*1e-6,lw=0.6,color='b',label='f 210 = ' + '%.2f' %(freq_modes[2][1][0]*1e-6) +' MHz' )
    #ax.axvline(x=freq_modes[2][1][1]*1e-6,lw=0.7,color='c',label='f 211 = ' + '%.2f' %(freq_modes[2][1][1]*1e-6) +' MHz' )
    #ax.axvline(x=freq_modes[2][2][0]*1e-6,lw=0.7,color='b',label='f 220 = ' + '%.2f' %(freq_modes[2][2][0]*1e-6) +' MHz' )
    #ax.axvline(x=freq_modes[2][2][1]*1e-6,lw=0.7,color='m',label='f 221 = ' + '%.2f' %(freq_modes[2][2][1]*1e-6) +' MHz' )
    #ax.axvline(x=freq_modes[2][2][2]*1e-6,lw=0.7,color='k',label='f 222 = ' + '%.2f' %(freq_modes[2][2][2]*1e-6) +' MHz' )
    ax.legend()
    fig.suptitle(title)
        

    frequency1=freq[0][peaks[0]]
    error1=abs(((frequency1-freq_modes[1][1][0])/(freq_modes[1][1][0])))
    print(error1)
    
    frequency2=freq[0][peaks[1]]
    error2=abs(((frequency2-freq_modes[1][1][1])/(freq_modes[1][1][1])))
    print(error2)
    
    frequency3=freq[0][peaks[2]]
    error3=abs(((frequency3-freq_modes[2][1][0])/(freq_modes[2][1][0])))
    print(error3)
    

    
##############################################################################
def Plot_Electric_Field(): 
    """The function Plot_Electric_Field() plots the summed total z-component
    transient electric field in the node volume boundend by the nodes labeld
    as:(i1_fld, j1_fld, k1_fld) to (i2_fld, j2_fld, k2_fld) which are 
    specified in the input parameter txt file."""
    
    title=('Transient Ez Electric Field, Node (16,16,16)-(16,16,17)')
    y_max=1.1*np.amax(Ez_out[0])   
    fig, ax =plt.subplots()
    ax.plot((time[0])*(10**9),Ez_out[0],label='Ez',color='k',lw=0.5)
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_ylim(-y_max, y_max)
    ax.set_xlim(0,simtime*(10**9))
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Electric Field (V/m)')
    fig.suptitle(title)
    
    return None
##############################################################################
def post_processing():
    """The function post_processing() is a function that calls on subfunctions
    to produce various plots and/or videos of the date outputed by the FDTD
    simulation performed by Time_advance()fucniton. The function 
    post_processing() is called at the end of the main() function when all 
    FDTD simulation functions have been performed"""
    
    '''Plot the transient Electric Field'''
    Plot_Electric_Field() 
    
    '''Compute and plot the fourer spectrum'''
    Fourier_Specturm();
    
    return None
##############################################################################
############################################################################## 
def main():
    
    '''Readin data for input Program Control txt file'''
    readin_parameter_data_file()
    
    '''Setup the simulation: Initiliaze Data Arrays and Constants'''
    setup_FDTD_sim()
    
    '''Perform time advacning FDTD update equations for Maxwells equations'''
    TimeAdvance()

    '''Perform post_proccessing of data to produce plots and video'''
    post_processing()
 
    return None
##############################################################################
##############################################################################
main()
##############################################################################
##############################################################################
##############################################################################
##############################################################################
