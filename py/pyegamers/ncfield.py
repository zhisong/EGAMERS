import numpy as np
import matplotlib.pyplot as plt

def read_ncfield(self, ncfield_name):
    fo = self.OutputData()

    from netCDF4 import Dataset
    import keyword

    rootgrp = Dataset(ncfield_name, "r", format="NETCDF4")

    for name in rootgrp.variables:
        if name in keyword.kwlist:
            name1 = name + '_'
        else:
            name1 = name
        setattr(fo, name1, np.array(rootgrp.variables[name][:]))
    
    re = rootgrp["/radial_element"]
    radial_attributes = ['lfieldoutput', 'igrid_type', 'XR1', 'XR2', 'SIG1', 'SIG2']
    for p in radial_attributes:
        setattr(fo, p, getattr(re, p))

    rootgrp.close()

    nele = fo.radial_element.size
    H1_elements = np.arange(1, nele-1, 2)
    H2_elements = np.concatenate([np.arange(0,nele,2), [nele-1]])
    fo.r = fo.radial_element[H2_elements]
    
    nax0 = np.zeros([fo.time.size, 1])

    fo.eta_H1 = np.concatenate([nax0, fo.eta[:, H1_elements], nax0], 1)
    fo.eta_H2 = fo.eta[:, H2_elements]

    fo.lambda_H1 = np.concatenate([nax0, fo.lambda_[:, H1_elements], nax0], 1)
    fo.lambda_H2 = fo.lambda_[:, H2_elements]

    return fo

def compute_ncfield(self, r=np.linspace(0,1,101), tslice=[0], ncfield=None):
    from scipy.interpolate import CubicHermiteSpline
    if ncfield==None:
        ncfield = self.ncfield

    tslice = np.atleast_1d(tslice)
    t = ncfield.time[tslice]

    lambda_CH = CubicHermiteSpline(ncfield.r, ncfield.lambda_H1[tslice,:], ncfield.lambda_H2[tslice,:], axis=1)
    eta_CH = CubicHermiteSpline(ncfield.r, ncfield.eta_H1[tslice,:], ncfield.eta_H2[tslice,:], axis=1)

    return t, r, lambda_CH(r), eta_CH(r)

def plot_field_radial(self, r=np.linspace(0,1,101), tslice=[0], variable='lambda', ncfield=None):

    tslice = np.atleast_1d(tslice)
    t, r, lambda_, eta = self.compute_ncfield(r=r, tslice=tslice, ncfield=ncfield)
    nt = tslice.size

    for i in range(nt):
        if variable == 'lambda':
            plt.plot(r, lambda_[i,:])
        else:
            plt.plot(r, eta[i,:])
        
    plt.legend(['t=' + "{0:.4g}".format(t1) for t1 in t ])

def plot_field_time(self, r=[0.5], tslice_start=0, tslice_end='last', variable='lambda', ncfield=None):
    r = np.atleast_1d(r)
    nr = r.size
    if ncfield==None:
        ncfield = self.ncfield

    if tslice_end == 'last':
        tslice_end = ncfield.time.size - 1

    t, r, lambda_, eta = self.compute_ncfield(r=r, tslice=np.arange(tslice_start, tslice_end), ncfield=ncfield)

    for i in range(nr):
        if variable == 'lambda':
            plt.plot(t, lambda_[:,i])
        else:
            plt.plot(t, eta[:,i])
        
    plt.legend(['r=' + "{0:.4g}".format(r1) for r1 in r ])
    

def plot_spectrogram(self, r=0.5, tslice_start=0, tslice_end='last', variable='lambda', ncfield=None, **kwargs):

    r = np.atleast_1d(r)

    if ncfield==None:
        ncfield = self.ncfield

    if tslice_end == 'last':
        tslice_end = ncfield.time.size - 1

    t, r, lambda_, eta = self.compute_ncfield(r=r, tslice=np.arange(tslice_start, tslice_end), ncfield=ncfield)

    from scipy import signal

    fs = (t.size - 1) / (t[-1] - t[0]) 

    if variable == 'lambda':
        data = lambda_[:,0]
    else:
        data =  eta[:,0]

    plt.specgram(data, Fs=fs, **kwargs)