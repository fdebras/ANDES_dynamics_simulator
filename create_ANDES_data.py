import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import interpolate
import pickle
import random
from parameters_and_functions_WASP76 import *


###################### Read the data and create interpolations
wlen_orders = np.loadtxt(dire+orders_file)  #SPIRou order wavelength

## Stellar data from phoenix
star_spec_read = fits.open(dire+file_star_spec) #stellar spectra
star_wl_read = fits.open(dire+file_star_wl) # wavelength of stellar spectra
star_spec = star_spec_read[0].data 
star_wl = star_wl_read[0].data

### Stellar flux interpolation
star_interp = interpolate.CubicSpline(star_wl, star_spec)



if phase_dependency:
    model_wl = np.loadtxt(dire_model_phase+prefix_model+str(phase_model[0])+suffix_model)[:,0]*1e9
    model_tdepth = []
    model_interp = []
    for i in range(len(phase_model)):
        model_tdepth.append(np.loadtxt(dire_model_phase+prefix_model+str(phase_model[i])+suffix_model)[:,1])
        if phase_model[i]>180:
            phase_model[i]-=360
    phase_model = np.array(phase_model)
    model_tdepth = np.array(model_tdepth)/100
    model_phase_interp = interpolate.CubicSpline(phase_model, model_tdepth)
    
else:
## Much simpler without the time dimension
    model_radius = np.loadtxt(dire+file_model_radius) # Planetary radius
    model_tdepth = (model_radius/(Rs*1e5))**2  # Transit depth
    model_wl = np.loadtxt(dire+file_model_wl) # Wavelentgh of planetary model
    
    ### transit depth interpolation
    model_interp = interpolate.CubicSpline(model_wl, model_tdepth)

### tellurics
tellurics_tot = np.loadtxt(dire+tellurics_file)
tellurics_wl = tellurics_tot[:,0]
tellurics_transm = tellurics_tot[:,1]

### tellurics interpolation
tellurics_interp = interpolate.CubicSpline(tellurics_wl, tellurics_transm)

### Read the observed HD 189 data
with open(dire+file_planet,'rb') as specfile:
    orders_spirou,wl_spirou,I_spirou,blaze_spirou,tellurics_spirou,T_obs_spirou,phase_spirou,\
        window_spirou,berv_spirou,vstar_spirou,airmass_spirou,SN_spirou = pickle.load(specfile)

## interpolation of the airmass for the tellurics
airmass_interp = interpolate.CubicSpline(T_obs_spirou,airmass_spirou)

#############################################################################################

######################## Prepare observations
### Given the SN scaling between SPIRou and ANDES, we increase the number of observations
### Is the scaling is 10, hence a 10 times higher SNR, we expect the same number of obs
### If the scaling is 1<n<10, we will have 10/n**2 more observation to ensure such a SNR with ANDES
T_obs = np.linspace(min(T_obs_spirou),max(T_obs_spirou),len(T_obs_spirou)*int((10/SN_scaling)**2))

### Adapt berv and airmass to this scaling
berv = np.linspace(berv_spirou[0],berv_spirou[-1],len(T_obs))
airmass = airmass_interp(T_obs)

### Calculate phase
phase  = (T_obs - T0)/Porb
phase -= int(phase[-1])  

### Compute transit window
flux     = compute_transit(Rp,Rs,ip,T0,ap,Porb,ep,wp,ld_mod,ld_coef,T_obs)
window       = (1-flux)/np.max(1-flux)

### Compute Planet-induced RV
Vs           = get_rvs(T_obs,Ks,Porb,T0)
Vc           = V0 + Vs - berv  #Geocentric-to-barycentric correction

#Compute planet RV
Vp = rvp(phase,Kp,V_inj)

### If we consider phase dependency in the model, create the array of interpolation
if phase_dependency:
    model_interp = []
    for i in range(len(T_obs)):
        if (phase[i]*360)<min(phase_model):
            model_phase = model_tdepth[0]
        elif (phase[i]*360)>max(phase_model):
            model_phase = model_tdepth[-1]
        else:
            model_phase = model_phase_interp(phase[i]*360)
        model_interp.append(interpolate.CubicSpline(model_wl, model_phase))
        



### Now we create the stellar data, looping over each SPIRou order
wl = []
data = []
SN = []
blaze = []
tellurics = []
orders_final = []
for i in range(len(orders_spirou)):
# for i in [25]:
    
    try:
        no=np.where(wlen_orders[:,0]==orders_spirou[i])[0][0]
    except:
        continue
    lmin = wlen_orders[no,1]
    lmax = wlen_orders[no,2]
    
    ### First, we prepare over resolved data that will be broadened afterwards
    ### We take slightly larger wavelength interval to ensure a correct interpolation
    ### for the broadening later on
    wl_order_high = [lmin*0.99]
    wl_next = wl_order_high[0]
    while wl_next <= lmax*1.01:
        wl_next += wl_next/(resolution*5)
        wl_order_high.append(wl_next)
    wl_order_high = np.array(wl_order_high)
    data_order_high = np.zeros((len(T_obs),len(wl_order_high)))
    data_order_high_tell = np.zeros((len(T_obs),len(wl_order_high)))
    
    ### calcualte oversampled tellurics
    if consider_tellurics:
        tell_high = tellurics_interp(wl_order_high).clip(min=0) ### ensure that there are no negative values
        
        

    
    ### Then we create the final wavelength range for our data, at the given resolution
    ### We do as in SPIRou : the pixel size is half the resolution
    wl_order = [lmin]
    wl_next = lmin
    while wl_next <= lmax:
        wl_next += wl_next/(2*resolution)
        wl_order.append(wl_next)
    wl_order = np.array(wl_order)
    
    ### Calculate an estimate of  final tellurics
    if consider_tellurics:
        tell_fin = tellurics_interp(wl_order).clip(min=0)
    
    data_order = np.zeros((len(T_obs),len(wl_order)))
    SN_order = np.zeros((len(T_obs)))

    ### We prepare the noise map at the given SNR
    noise_order = np.zeros(data_order.shape)
    for t in range(len(data_order)):
        noise_order[t] = np.random.normal(scale=1./np.mean(SN_spirou[i])/SN_scaling,size=np.shape(data_order)[1])
        
    ### IF we want blaze, we interpolate the SPIRou blaze and add it to the noise
    if consider_blaze:
        blaze_interp = interpolate.CubicSpline(wl_spirou[i], blaze_spirou[i][0])
        blaze_order = blaze_interp(wl_order)
        
        blaze_and_noise_order = 1+(noise_order)*np.max(np.sqrt(blaze_order))/np.sqrt(blaze_order)
    

    ### We start by interpolating the star considering berv and planet induced RV
    for t in range(len(T_obs)):
        data_order_high[t] = star_interp(wl_order_high*10*(1-Vc[t]/c0)) ### *10 because phoenis is in Angstrom
        data_order_high[t] /= np.max(data_order_high[t]) ### We normalize the data to 1, not a big deal
        
        ### In SPIRou, the flux can vary a lot from one observation to another because of imperfect
        ### centering of the fiber. I don't know how does that trnaslate to the ELT. You can leave or suppress that
        if consider_flux_variations:
            data_order_high[t]*= random.gauss(np.mean(I_spirou[i]),np.std(np.mean(I_spirou[i],axis=1)))
            
        ### Add the planet, doppler shifted by its orbital velocity
        if window[t]>0:
            if phase_dependency:
                data_order_high[t] *= (1.-model_interp[t](wl_order_high*(1-(Vc[t]+Vp[t])/c0)))
            else:
                data_order_high[t] *= (1.-model_interp(wl_order_high*(1-(Vc[t]+Vp[t])/c0))*window[t])
            
        ### Add tellurics
        if consider_tellurics:
            data_order_high_tell[t] = data_order_high[t]*tell_high**airmass[t]
            
            ### Broadening of the new data at the instrumental resolution
            wl_order_high_broad,data_order_high_broad = broaden(wl_order_high,data_order_high_tell[t],sigma*1000)
            ### Interpolation of this broadened model
            data_high_broad_interp = interpolate.CubicSpline(wl_order_high_broad, data_order_high_broad)
        else:
            wl_order_high_broad,data_order_high_broad = broaden(wl_order_high,data_order_high[t],sigma*1000)
            ### Interpolation of this broadened model
            data_high_broad_interp = interpolate.CubicSpline(wl_order_high_broad, data_order_high_broad)
        
        ### Integrate the intensity over the pixel size. This is performed very naively here
        tmp = np.zeros(len(wl_order))
        for pp in pixel_array: 
            tmp += data_high_broad_interp(wl_order/(1.0+(pp)/c0))
        data_order[t] = tmp/len(pixel_array)
        
    ### Add the noise, eventually including blaze, and calculate SNR at the center of the order
    mid = int(len(wl_order)/2)
    if consider_blaze:
        data_order*= blaze_and_noise_order
        SN_order =1./ np.std(blaze_and_noise_order[:,mid-int(mid/4):mid+int(mid/4)],axis=1)
    else:
        data_order *= 1+noise_order
        SN_order = 1./np.std(noise_order[:,mid-int(mid/4):mid+int(mid/4)],axis=1)


    ### Store the data
    orders_final.append(orders_spirou[i])
    data.append(data_order)
    wl.append(wl_order)
    SN.append(SN_order)
    blaze.append(np.tile(blaze_order,(len(T_obs),1))) 
    ### To be consistent with SPIRou fits file we do this but it can be changed
    tellurics.append((np.tile(tell_fin,(len(T_obs),1)).T**airmass).T)
    

### Save the data
savedata = (orders_final,wl,data,blaze,tellurics,T_obs,phase,window,berv,V0+Vs,airmass,SN)
with open(dire+save_file, 'wb') as specfile:
    pickle.dump(savedata,specfile)
print("DONE")

