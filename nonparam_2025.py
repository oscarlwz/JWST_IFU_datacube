import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit as curvefit
from  scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy import interpolate
from scipy.interpolate import CubicSpline

from scipy.special import ndtr as ndtr
import pdb
import os
from astropy.cosmology.funcs import *
from matplotlib import pyplot as plt
from copy import copy as copy

import sys
#sys.path.append('/Users/weizheliu/cos-ulirg/')
import pickle

#==========================================================================


def multigauss(p, x, data, err):
    model = 0
    pp = p
    keyp = list(pp.keys())
    ng = int(np.size(keyp)/3)
    #pdb.set_trace()
    for i in range(ng):
        model = model+gauss(x,p[keyp[3*i]],p[keyp[3*i+1]],p[keyp[3*i+2]])
    return (data - model)/err



c = 299792.458

def wav2vel(wav,wav0):
    cspeed = 2.99792458e5
    z = (wav/wav0-1)
    vel = ((1+z)**2-1) / ((1+z)**2+1)
    return vel*cspeed
def vel2wav(vel,wav0):
    cspeed = 2.99792458e5
    b = vel/cspeed
    wav = wav0*((1+b)/(1-b))**0.5
    return wav

def z2vel(z):
    cspeed = 2.99792458e5
    vel = ((1+z)**2-1) / ((1+z)**2+1)
    return vel*cspeed


def vel2z(vel):
    cspeed = 2.99792458e5
    b = vel/cspeed
    z = ((1+b)/(1-b))**0.5-1
    return z


def fwav2fvel(wav,wav0):
    cspeed = 2.99792458e5
    z = (wav/wav0-1)
    fscale = 4*(1+z)/ (2+2*z+z**2)**2
    fscale = wav0/(fscale*cspeed)
    return fscale  #can be multiplied to f_lambda



def fitfun(x,p0,p1,p2,p3,p4):
    return p0+p1/(pow(x-p3,p2)+p4)

def powfit(x,alpha=1.0,scale=1.0):
    return scale*pow(x,alpha)

def linefit(x,alpha=1.0, offset=-3):
    return alpha*x + offset

def gaussint(x):
    integ = np.zeros(np.size(x))
    err = np.zeros(np.size(x))
    i = 0
    for x0 in x:
        #integ[i], err[i] = integrate.quad(gauss,-1.*np.inf,x0,args=(0.,1.,1./((2*np.pi)**0.5)))
        #i = i +1
        integ[i] = ndtr(x0)
        i = i+1
        #print('haha',x0,integ[i])
    return integ


#mpc = 3.086e24
#
def luminosity(flux,err,z,log=False):
    from astropy.cosmology import WMAP9 as cosmo
    mpc = 3.086e24
    dis = cosmo.luminosity_distance(z).value*mpc
    if log:
        return flux+np.log10(4*np.pi*dis**2), err+np.log10(4*np.pi*dis**2)
        #return flux+np.log10(4*np.pi*dis**2), np.log10(10**(flux+err)-10**flux)+np.log10(4*np.pi*dis**2)
    else:
        return flux*4*np.pi*dis**2,err*4*np.pi*dis**2




#all the wave, flux are assumed to be in the rest-frame
class NonParam:
    def __init__(self,wave,flux,error, wave0=1216., wav_nonpar=[], w_i=-999., w_j=-999., ew_ind=[], useraw=True, eps=0.001,contipar=[],contipar_err=[],renorm=None, contitype='powerlaw', n_continuum=2,continuumflag='line',wav_chazhi=[],waicha=[],absorption=False,con_err=[]):
        self.con_err = con_err
        self.useraw = useraw
        #pdb.set_trace()
        #if np.size(flux) > 5:
        #    self.valid = True
        #else:
        #    self.valid = False
        #    retrun
        if np.size(wav_nonpar) > 0:
            w_i = wav_nonpar[0]
            w_j = wav_nonpar[1]
            ind = np.where((wave > w_i) & (wave < w_j))    
        elif w_i > 0.:
            ind = np.where((wave > w_i) & (wave < w_j))
        else:
            ind = np.arange(np.size(wave))

        self.valid = True
        if np.size(ind) < 1:
            self.valid =False
            #pdb.set_trace()
            #return None
        #print('indddda',ind)
        #if np.size(ew_ind) > 1: 
        #    ind = ind_lya 
        flux_chazhi_all = []
        wave_chazhi_all = []
        ncha = []
        flux_origin = copy(flux)
        if np.size(wav_chazhi)>0:
            for i in np.arange(np.size(wav_chazhi)/4):
                #for w_cha in wav_chazhi:
                w_cha = wav_chazhi
                chazhi1 = np.where((wave > w_cha[int(0+i*4)]) & (wave < w_cha[int(1+i*4)]))
                chazhi2 = np.where((wave > w_cha[int(2+i*4)]) & (wave < w_cha[int(3+i*4)]))
                wavnew = np.append(wave[chazhi1],wave[chazhi2])
                fluxnew = np.append(flux[chazhi1],flux[chazhi2])
                p = np.polyfit(wavnew,fluxnew,2)
                # cannot be used in the current scipy, again ...
                #f = interpolate.interp1d(np.append(wavnew,fluxnew,kind='cubic'))
                chazhi = np.where((wave > w_cha[int(1+i*4)]) & (wave < w_cha[int(2+i*4)]))
                flux[chazhi] = np.polyval(p,wave[chazhi])
                flux_chazhi_all = np.append(flux[chazhi],flux_chazhi_all)
                wave_chazhi_all = np.append(wave[chazhi],wave_chazhi_all)
                ncha = np.append(np.size(chazhi),ncha)
        wave_waicha_all = []
        flux_waicha_all = []
        if np.size(waicha) > 0:
            for i in np.arange(np.size(waicha)/3):  
                waichazhi1 = np.where((wave > waicha[int(0+i*3)]) & (wave < waicha[int(1+i*3)]))
                waichazhi2 = np.where((wave > waicha[int(1+i*3)]) & (wave < waicha[int(2+i*3)]))
                waichazhiall = np.where((wave > waicha[int(0+i*3)]) & (wave < waicha[int(2+i*3)]))
                fmax = np.max(flux[waichazhi1])
                p2,p2cov = curvefit(gauss,wave[waichazhi1],flux[waichazhi1]-flux[waichazhi1][0],p0=[1215.6,10,fmax])
                offset = flux[waichazhi2][0] - gauss(wave[waichazhi2],p2[0],p2[1],p2[2])[0] - flux[waichazhi1][0]
                flux_waicha_all = np.append(flux[waichazhi2],flux_waicha_all)
                wave_waicha_all = np.append(wave[waichazhi2],wave_waicha_all)
        
        self.wave_input = wave[ind]
        self.flux_input = flux[ind]
        self.wave_chazhi_all = wave_chazhi_all
        self.flux_chazhi_all = flux_chazhi_all
        self.ncha = ncha
        self.flux_waicha_all = flux_waicha_all
        self.wave_waicha_all = wave_waicha_all
        self.error0 = error
        self.wave = wave[ind]
        self.vel = wav2vel(wave[ind],wave0) 
        self.wave0 = wave0  
        self.flux = flux[ind]
        self.error = error[ind]
        self.contitype=contitype
        self.contipar = contipar
        self.contipar_err = contipar_err
        self.renorm = renorm
        self.absorption = absorption
    
    #eps is chosen relatively large now, given the step size is not very small
    def v(self, vthres = 0.2, eps=0.1, plot=True, vkms=None,uplim=False,dowtavg=False,weq_uplim=False):    
        con_err = copy(self.con_err)
        #print('wolegequ  ',self.contipar)
        #breakpoint()
        if self.useraw:
            vel = self.vel
            #pdb.set_trace()
            wave = self.wave
            #renorm = np.median(self.flux_input)
            #renorm = copy(self.renorm)
            flux = self.flux_input# / renorm
            err = self.error# / renorm
            contipar = copy(self.contipar)
            contipar_err = copy(self.contipar_err)
            if self.contitype == 'powerlaw':
                try:
                    contipar[1] = contipar[1] #/renorm
                    contipar_err[1] = contipar_err[1] #/renorm #normalized
                    continuum = contipar[1]*pow(wave,contipar[0])
                except:
                    pdb.set_trace()
                self.continuum = continuum  #*renorm    
            elif self.contitype == 'spline':
                continuum = contipar[0](wave)
                self.continuum = continuum
            elif self.contitype == 'poly2' or self.contitype == 'poly':
                continuum = np.polyval(contipar,wave)
                self.continuum = continuum
            else:
                pdb.set_trace()

            #if self.contitype != 'spline' and  self.contitype != 'powerlaw':
            #    pdb.set_trace()
            flux_nocon = flux - continuum
            if self.absorption:
                flux_nocon = continuum - flux
            #print(np.size(continuum),np.size(flux))
            #print(np.size(flux_nocon),np.size(wave))
            if vkms == None:
                # Compute cumulative flux
                cumulative_flux = np.cumsum(flux_nocon)
                total_flux = cumulative_flux[-1]

                # Calculate the 10% and 90% cumulative flux levels
                flux_vthres = vthres * total_flux

                # Interpolate to find the velocities at 10% and 90% flux levels
                v_out0 = np.interp(flux_vthres, cumulative_flux, vel)



                flux_t = integrate.simps(flux_nocon,wave)
                n_wave = np.size(wave)
                dv1 = vel[1] - vel[0]
                flux_tmp_err2 = 0. #flux error squared
                flux_tmp_con_err2 = 0.
                for i in np.arange(n_wave-2)+1:
                    flux_tmp = integrate.simps(flux_nocon[:i],wave[:i])
                    #add up the errors for each flux bin (d_wave*err)
                    #err2 means it is error squared
                    flux_tmp_err2 = flux_tmp_err2 + (err[i]*(wave[i+1] - wave[i]))**2 
                    #errors consider the uncertainties in continuum
                    # ignore these uncertainties for now
                    if self.contitype == 'powerlaw':
                        con_err2 = (pow(wave[i],contipar[0])*contipar_err[1])**2 + (contipar[0]*contipar[1]*pow(wave[i],contipar[0]-1.)*contipar_err[0])**2
                    elif self.contitype == 'poly2':
                        con_err2 = wave[i]**4*contipar_err[0]**2+wave[i]**2*contipar_err[1]**2+contipar_err[2]**2
                    elif self.contitype == 'poly':
                        con_err2 = wave[i]**2*contipar_err[0]**2+contipar_err[1]**2
                    else:
                        con_err2 = 0 # to be modified to real errors
                        
                    flux_tmp_con_err2 = flux_tmp_con_err2 + con_err2*(wave[i+1] - wave[i])**2
                    test = flux_tmp/flux_t - vthres
                    


                    #print(flux_tmp, test,' !!!!?????')
                    if abs(test) < eps:
                        flux_tmp_2 = integrate.simps(flux_nocon[:i+1],wave[:i+1])
                        test2 = flux_tmp_2/flux_t - vthres
                        #print(test2,test,'dasdadasd')
                        if test2*test < 0:
                            
                            if abs(test2) < abs(test):
                                v_out = vel[i+1]
                                flux_out = flux_tmp_2
                                flux_tmp_err2 = flux_tmp_err2+(err[i+1]*(wave[i+2] - wave[i+1]))**2
                                if self.contitype == 'powerlaw':
                                    con_err2 = (pow(wave[i+1],contipar[0])*contipar_err[1])**2 + (contipar[0]*contipar[1]*pow(wave[i+1],contipar[0]-1.)*contipar_err[0])**2
                                else:
                                    con_err2 = 0
                                flux_tmp_con_err2 = flux_tmp_con_err2 + con_err2*(wave[i+2] - wave[i+1])**2
                                #print(flux_tmp_err2,' ?????')
                                j = i+2
                                j_con = i+2
                            else:
                                v_out = vel[i]
                                flux_out = flux_tmp
                                flux_tmp_err2 = flux_tmp_err2+(err[i]*(wave[i+1] - wave[i]))**2
                                if self.contitype == 'powerlaw':
                                    con_err2 = (pow(wave[i],contipar[0])*contipar_err[1])**2 + (contipar[0]*contipar[1]*pow(wave[i],contipar[0]-1.)*contipar_err[0])**2
                                else:
                                    con_err2 = 0
                                flux_tmp_con_err2 = flux_tmp_con_err2 + con_err2*(wave[i+1] - wave[i])**2
                                #print(flux_tmp_err2,' !!!!')
                                j = i+1
                                j_con = i+1
                            flux_tmp_err = flux_tmp_err2**0.5
                            flux_tmp_con_err = flux_tmp_con_err2**0.5
                            #The error in v_out is the velocity difference between flux_tmp & flux_tmp+flux_tmp_err
                            while(j < n_wave-1):
                                flux_tmp_uperror = integrate.simps(flux_nocon[:j],wave[:j])
                                #print(n_wave,j,flux_tmp_uperror,flux_tmp,flux_tmp_err)
                                if abs(flux_tmp_uperror-flux_tmp) > flux_tmp_err:
                                    vout_err = abs(vel[j] - v_out) 
                                    break
                                j = j+1
                                #print(j)
                            # or considering the error in the fitted continuum: vout_con_err
                            if j == n_wave-1:
                                vout_err = abs(vel[-1] - v_out)
                            while(j_con < n_wave-1):
                                flux_tmp_con_uperror = integrate.simps(flux_nocon[:j_con],wave[:j_con])     

                                if abs(flux_tmp_con_uperror-flux_tmp) > flux_tmp_con_err:
                                    vout_con_err = abs(vel[j_con] - v_out)
                                    #if vthres == 0.5:
                                    #    pdb.set_trace()
                                    break
                                j_con = j_con+1
                            #the maximum error is the end velocity - vout, when error is too large
                            if j_con == n_wave-1:
                                vout_con_err = abs(vel[-1] - v_out)
                            #vout_err = abs(vel[i+1]-vel[i])
                            #why flux_err and vout_err so small sometimes
                            if vout_err <= 0:
                                vout_err = vel[j]-vel[j-1]
                            if vout_con_err <= 0:
                                vout_con_err = vel[j]-vel[j-1]
                            #print('!!!!!!')
                            #return v_out, vout_err, vout_con_err, flux_out, flux_tmp_err, flux_tmp_con_err
                            #plot the regions for the 
                            #return v_out, vout_err, vout_con_err, flux_out*renorm, flux_tmp_err*renorm, flux_tmp_con_err*renorm    
                            return -996, vout_err, vout_con_err, flux_out, flux_tmp_err, flux_tmp_con_err
            else:
                #if vkms0 == None:
                #    vkms0 == -20000
                ind = np.where(vel < vkms)# and vel > vkms0)
                if np.size(ind) > 1:
                    fluxkms = integrate.simps(flux_nocon[ind],wave[ind])
                    errkms2 = np.sum(err[ind]**2)
                    errkms = errkms2**0.5
                else:
                    fluxkms = np.nan
                    errkms2 = np.nan
                    errkms = np.nan
                return fluxkms,errkms
                #    self.vnp = self.vel[i]

        else:
            vel = self.vel
            #pdb.set_trace()
            con_err = copy(self.con_err)
            if np.size(con_err) == 0:
                con_err = np.zeros(np.size(self.wave))
            wave = self.wave
            flux = self.flux_input 
            err = self.error
            flux_nocon = flux.copy()
            flux_nocon_err = err.copy()
            flux_nocon_vel = flux_nocon * fwav2fvel(wave,self.wave0)
            flux_nocon_vel_err = flux_nocon_err * fwav2fvel(wave,self.wave0)
            minid = np.argmin(np.abs(flux_nocon-1.))
            flux_nocon_vel = flux_nocon_vel / fwav2fvel(wave[minid],self.wave0)
            flux_nocon_vel_err = flux_nocon_vel_err / fwav2fvel(wave[minid],self.wave0)
            if self.absorption:
                flux_nocon = np.zeros(np.size(flux_nocon))+1 - flux_nocon
                flux_nocon_vel = np.zeros(np.size(flux_nocon_vel))+1 - flux_nocon_vel   
                if dowtavg:
                    if weq_uplim:
                        linewindow = np.where(np.abs(wave-self.wave0) < 0.5)
                        ew_lt0 = np.std(flux_nocon[linewindow])
                        #gamma = 
                        #wav0 =
                        b = 50*2**0.5
                        tp = 5
                        #w1 = (2*b/cspeed)**2 * np.log(tau0/np.log(2))
                        #w2 = b/c*gamma*wav0/c*(tau0-1.25393)/(np.pi**2)  #for tau > 1.254, Draine's book 9.27, P81
                        # use this formula for simplicity for tau < 1.254
                        cspeed = 2.99792458e5
                        w1 = np.pi**0.5*b/cspeed
                        w2 = tp/(1+tp/2/2**0.5)
                        ew_lt0_A = self.wave0 *(w1*w2) * ew_lt0
                        #ew_lt0 = np.nan
                        #gauss_line(50,0,nu,nu0)
                        #1-cf*(1-np.exp(-tau.astype(np.float)))
                        v_wtavg = np.nan
                        s_wtavg = np.nan
                        ew_lt0_err = np.nan
                        ew_lt0_A_err =np.nan
                        v_wtavg_err = np.nan
                        s_wtavg_err = np.nan
                        #print(ew_lt0_A,v_wtavg,s_wtavg)
                        #breakpoint()
                        return ew_lt0_A,v_wtavg,s_wtavg,ew_lt0
                    else:
                        ew_lt0 = 0.
                        ew_lt0_A = 0.
                        v_wtavg = 0.
                        s_wtavg = 0.
                        ew_lt0_err2 = 0.
                        ew_lt0_A_err2 = 0.
                        v_wtavg_err2 = 0.
                        s_wtavg_err2 = 0.
                        for iv in range(np.size(vel)-1):
                            if vel[iv] < 1000 and vel[iv] > -30000:  # now the criterion is not forced to be below 0km/s
                            #if vel[iv] < 0 and vel[iv] > -30000:
                                ew_lt0 = ew_lt0 + flux_nocon_vel[iv] * (vel[iv+1]-vel[iv])
                                ew_lt0_A = ew_lt0_A + flux_nocon[iv] * (wave[iv+1]-wave[iv])
                                v_wtavg = v_wtavg + vel[iv]*flux_nocon_vel[iv] * (vel[iv+1]-vel[iv])
                                ew_lt0_err2 = ew_lt0_err2 + (flux_nocon_vel_err[iv] * (vel[iv+1]-vel[iv]))**2
                                ew_lt0_A_err2 = ew_lt0_A_err2 + (flux_nocon_err[iv] * (wave[iv+1]-wave[iv]))**2
                                v_wtavg_err2 = v_wtavg_err2 + (vel[iv]*(flux_nocon_vel_err[iv] * (vel[iv+1]-vel[iv])))**2
                        v_wtavg = v_wtavg/ew_lt0
                        ew_lt0_err = ew_lt0_err2**0.5
                        ew_lt0_A_err = ew_lt0_A_err2**0.5
                        v_wtavg_err = v_wtavg_err2**0.5/ew_lt0
                        for iv in range(np.size(vel)-1):
                            if vel[iv] < 0 and vel[iv] > -30000:
                                s_wtavg = s_wtavg + (vel[iv]-v_wtavg)**2*flux_nocon_vel[iv] * (vel[iv+1]-vel[iv])
                                s_wtavg_err2 = s_wtavg_err2 + ((vel[iv]-v_wtavg)**2 * flux_nocon_vel_err[iv] * (vel[iv+1]-vel[iv]))**2
                        s_wtavg = (np.abs(s_wtavg/ew_lt0))**0.5
                        #breakpoint()
                        s_wtavg_err = 0.5*(s_wtavg_err2**0.5)/((np.abs(s_wtavg*ew_lt0))**0.5)
                
            #print(np.size(continuum),np.size(flux))
            #print(np.size(flux_nocon),np.size(wave))
            if vkms == None:
                flux_t = integrate.simps(flux_nocon,wave)
                n_wave = np.size(wave)
                dv1 = vel[1] - vel[0]
                flux_tmp_err2 = 0. #flux error squared
                flux_tmp_con_err2 = 0.

                vout_samples = []
                noisy_flux_all = []
                flux_err = (err**2+con_err**2)**0.5
                num_iterations=10000
                
                # inspired by chatgpt ...  easier by my original code
                for _ in range(num_iterations):
                    # Add noise to the flux based on the flux error
                    noisy_flux = flux_nocon + np.random.normal(0, flux_err)
                    noisy_flux_all.append(noisy_flux)
                    # Compute cumulative flux

                    # Compute cumulative flux
                    cumulative_flux = np.cumsum(noisy_flux)
                    infinite_id = np.where(np.isfinite(cumulative_flux) == False)
                    if np.size(infinite_id) > 0: 
                        cumulative_flux = np.nancumsum(noisy_flux)
                    total_flux = cumulative_flux[-1]

                    # Calculate the 10% and 90% cumulative flux levels
                    flux_vthres = vthres * total_flux

                    # Interpolate to find the velocities at 10% and 90% flux levels
                    v_out_tmp = np.interp(flux_vthres, cumulative_flux, vel)
                    vout_samples.append(v_out_tmp)
                vout_0 = np.nanmean(vout_samples)
                vout_err0 = np.nanstd(vout_samples)
                vout_con_err, flux_out, flux_tmp_err, flux_tmp_con_err = 0, 0, 0, 0
                    
                return vout_0, vout_err0, vout_con_err, noisy_flux_all, flux_tmp_err, flux_tmp_con_err


            else:
                ind = np.where(vel < vkms)# and vel > vkms0)
                if np.size(ind) > 1:
                    fluxkms = integrate.simps(flux_nocon[ind],wave[ind])
                    errkms2 = np.sum(err[ind]**2)
                    errkms = errkms2**0.5
                else:
                    fluxkms = np.nan
                    errkms2 = np.nan
                    errkms = np.nan
                if dowtavg and self.absorption:
                    return ew_lt0_A,v_wtavg,s_wtavg,fluxkms,errkms
                else:
                    return fluxkms,errkms



    
    def plot(self,vlist=[],vnamelist=[],outpath='',targname='',normIC=False):
        self.outpath = outpath
        self.targname = targname
        vel = self.vel
        flux_vel = copy(self.flux_input)/c*copy(self.wave0)
        self.flux_vel = flux_vel
        continuum_vel = copy(self.continuum)/c*copy(self.wave0)
            
        fig = plt.figure(figsize=(8,6), dpi=150)
        ax = plt.subplot(111)
        #plt.margins(0.1)
        plt.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.96,
                            wspace=0.15, hspace=0.28)
        MAX_FLUX_all = 1.1*np.max(flux_vel)
        from plotcos_ulirg_bayesvp import plot_spectrum 
        plot_spectrum('', vel, flux_vel, -999, -999, 0, MAX_FLUX_all, plotwindow = False,  error=self.error, smooth=0, ax0 = ax, linemarker=0, linename = '', z=0., linenameout=False, plotvel=True)
        #plot_spectrum('', self.vel, self.flux, -999, -999, 0, MAX_FLUX_all, plotwindow = True, window=window, wc=wc, labeltext=labeltext, error=self.flux_e_rest, smooth=0, ax0 = ax1, linemarker=482, linename = 'all', z=self.z, linenameout=False)
        ax.plot(self.vel,continuum_vel,'g')
        co = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        for i in range(np.size(vlist)):
            ax.axvline(x=vlist[i],linestyle='--',color=co[i],label=vnamelist[i])
           
        fig.tight_layout()
        plt.legend(loc=0)
        #pdb.set_trace()
        plt.savefig(outpath+targname+'_vlist.png')
        plt.close(fig)
        if normIC:
            norm_id = np.argmin(np.abs(np.absself.wave-self.wave0))
            flux_normed = copy(flux_vel)/flux_vel[norm_id]
            max_normed = np.max(flux_normed)
            fig = plt.figure(figsize=(8,6), dpi=150)
            ax = plt.subplot(111)
            plt.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.96,
                            wspace=0.15, hspace=0.28)
            plot_spectrum('', vel, flux_normed, -999, -999, 0, max_normed*1.1, plotwindow = False,  error=self.error, smooth=0, ax0 = ax, linemarker=0, linename = '', z=0., linenameout=False, plotvel=True)
        fig.tight_layout()
        plt.legend(loc=0)
        plt.savefig(outpath+targname+'_normedvlist.png')
        plt.close(fig)


    def ew(self,ew_ind=[],wav_ew=[],con_window=None,useraw=True,contipar=None,contipar_err=None,z=-999.,uplim=False,con_type=''):
        con_err = copy(self.con_err)
        if np.size(ew_ind) > 0:
            ew_ind = np.where((self.wave_input > wav_ew[0]) & (self.wave_input < wav_ew[1])) 
            ew_wave = self.wave_input[ew_ind]
            ew_flux = self.flux_input[ew_ind]
            ew_fluxerr = self.error0[ew_ind]
        else:
            ew_wave = self.wave_input
            ew_flux = self.flux_input
            ew_fluxerr = self.error0
        if self.useraw:
            if uplim == True:
                dw = ew_wave[1] - ew_wave[0]
                lineflux = 5*(np.sum((ew_fluxerr*dw)**2))**0.5  # 5 sigma limit
                ew=ew_err=ew_con_err=lineflux_err=lineflux_con_err = np.nan
            else:
                #if np.size(contipar) == 2 :
                if self.contitype == 'powerlaw':
                    continuum = contipar[1]*pow(ew_wave,contipar[0])
                    con_err2 = (pow(ew_wave,contipar[0])*contipar_err[1])**2 + (contipar[0]*contipar[1]*pow(ew_wave,contipar[0]-1.)*contipar_err[0])**2
                elif self.contitype == 'spline':
                    continuum = contipar[0](ew_wave)
                    self.continuum = continuum
                elif self.contitype == 'poly2':
                    continuum = np.polyval(contipar,ew_wave)
                    con_err2 = ew_wave**4*contipar_err[0]**2+ew_wave**2*contipar_err[1]**2+contipar_err[2]**2
                    self.continuum = continuum
                elif self.contitype == 'poly':
                    continuum = np.polyval(contipar,ew_wave) 
                    con_err2 = ew_wave**2*contipar_err[0]**2+contipar_err[1]**2
                    self.continuum = continuum
                else:
                    pdb.set_trace()

                    con_err2 = 0  # to be modified
                con_err =con_err2**0.5
                ew = 0.
                ew_err2 = 0.
                ew_con_err2 = 0.
                lineflux = 0.
                lineflux_err2 = 0.
                lineflux_con_err2 = 0.
                n_wave = np.size(ew_wave)
                for i in np.arange(n_wave-1):
                    dw = ew_wave[i+1] - ew_wave[i]
                    ew = ew + dw * (ew_flux[i] - continuum[i])/continuum[i]
                    print('hahah',(ew_flux[i] - continuum[i])/continuum[i], continuum[i],dw,ew)
                    ew_con_err2 = ew_con_err2 + dw**2 * ((ew_fluxerr[i]/continuum[i])**2 + (ew_flux[i]*con_err[i]/(continuum[i]**2))**2)
                    ew_err2 = ew_err2 + (dw*ew_fluxerr[i]/continuum[i])**2
                    lineflux = lineflux + dw * (ew_flux[i] - continuum[i])
                    print('yi>>>>><<',dw,(ew_flux[i] - continuum[i]))
                    #lineflux_err2 = lineflux_err2 + (dw * ew_fluxerr[i])**2
                    #lineflux_con_err2 = lineflux_con_err2 + dw**2 * (ew_fluxerr[i]**2 + con_err[i]**2)
                    lineflux_con_err2 = lineflux_con_err2 + dw**2 * con_err[i]**2
                    lineflux_err2 = lineflux_err2 + (dw * ew_fluxerr[i])**2 - dw**2 * con_err[i]**2
                    #print(dw,ew,ew_err2,ew_flux[i],ew_fluxerr[i],continuum[i])
                    #pdb.set_trace()
                ew_err = ew_err2 ** 0.5
                lineflux_err = lineflux_err2 ** 0.5
                ew_con_err = ew_con_err2 ** 0.5
                lineflux_con_err = lineflux_con_err2 ** 0.5
                #breakpoint()
            return ew,ew_err,ew_con_err,lineflux,lineflux_err,lineflux_con_err

        # use continuum normalized data
        else:
            if uplim == True:
                dw = ew_wave[1] - ew_wave[0]
                lineflux = 5*(np.sum((ew_fluxerr*dw)**2))**0.5  # 5 sigma limit
                ew=ew_err=ew_con_err=lineflux_err=lineflux_con_err = np.nan
            else:
                if np.size(con_err) == 0:
                    con_err = np.zeros(np.size(ew_flux))
                #if np.size(contipar) == 2 :
                ew = 0.
                ew_err2 = 0.
                ew_con_err2 = 0.
                lineflux = 0.
                lineflux_err2 = 0.
                lineflux_con_err2 = 0.
                n_wave = np.size(ew_wave)
                for i in np.arange(n_wave-1):
                    dw = ew_wave[i+1] - ew_wave[i]
                    ew = ew + dw * (ew_flux[i])
                    ew_con_err2 = ew_con_err2 + dw**2 * ((ew_fluxerr[i])**2 + (ew_flux[i]*con_err[i])**2)
                    ew_err2 = ew_err2 + (dw*ew_fluxerr[i])**2
                    lineflux = lineflux + dw * (ew_flux[i])
                    #lineflux_err2 = lineflux_err2 + (dw * ew_fluxerr[i])**2
                    #lineflux_con_err2 = lineflux_con_err2 + dw**2 * (ew_fluxerr[i]**2 + con_err[i]**2)
                    lineflux_con_err2 = lineflux_con_err2 + dw**2 * con_err[i]**2
                    lineflux_err2 = lineflux_err2 + (dw * ew_fluxerr[i])**2 - dw**2 * con_err[i]**2
                    #if z>0.:
                    #    ew = ew/(1+z)
                    #    ew = ew_err/(1+z)
                ew_err = ew_err2 ** 0.5
                lineflux_err = lineflux_err2 ** 0.5
                ew_con_err = ew_con_err2 ** 0.5
                lineflux_con_err = lineflux_con_err2 ** 0.5
            #calclate the upper limits 
            return ew,ew_err,ew_con_err,lineflux,lineflux_err,lineflux_con_err


if __name__ == '__main__':

    ppxf_example_kinematics_sdss()
    import matplotlib.pyplot as plt
    plt.pause(1)


