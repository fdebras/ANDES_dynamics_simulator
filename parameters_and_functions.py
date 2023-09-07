
import numpy as np
import batman
from scipy import interpolate
from scipy import signal

################################### Orbital parameters ##############

### Ephemerides (to compute orbital phase)
T0       =  2458334.990899        #Mid-transit time [BJD]
Porb     = 2.218577                      #Orbital period [d]

### Transit parameters -- Compute the transit window
### Using batman python package https://lweb.cfa.harvard.edu/~lkreidberg/batman/
### Get the limb-darkening coefficients in H band from Claret+2011: https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/529/A75
Rp       = 79838.362 #Planet radius  [km]
Rs       = 536181.8  #Stellar radius [km] 
ip       = 85.0    #Transit incl.  [deg]
ap       = 8.67   #Semi-maj axis  [R_star]
ep       = 0.0     #Eccentricity of Pl. orbit
wp       = 0.0     #Arg of periaps [deg]
ld_mod   = "linear"     #Limb-darkening model ["nonlinear", "quadratic", "linear"]
ld_coef  = [0.35]  #Limb-darkening coefficients 


### Radial velocity info
Ks        = 0.2    # RV semi-amplitude of the star orbital motion due to planet [km/s]
Kp        = 151.4  # RV semi-aplitude of the planet
V_inj     = -5.0   # injected Doppler shift of the model
V0        = -2.2   #S tellar systemic velocity [km/s]
#####################################################################



############ Spectral resolution parameters #######
resolution = 115000
c0 = 3e5

sigma = c0/resolution/(2*np.sqrt(2*np.log(2)))

pixel_size = c0/resolution/2
npixel = 10
pixel_array     = np.linspace(-pixel_size/2,pixel_size/2,npixel) # Window for integrating over a SPIRou pixel [km/s] -- Velocity bins of 2.28 km/s (Donati+2020b)
##################################################

############# Signal to noise parameters #########
SN_scaling = 10
consider_flux_variations = True
consider_blaze = True
##################################################

############# Consider tellurics ? #########
consider_tellurics = True
##################################################


###################################### File names #####################
dire = "/home/adminloc/Bureau/Atmospheres/ELT/"
file_star_spec = "lte05000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
file_star_wl = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
file_planet = "HD189_SPIRou.pkl"

file_model_radius = "Rp_HD189_onlyH2O-exomol-VMR3-T900.txt"
file_model_wl = "lambdas_HD189_onlyH2O-exomol-VMR3-T900.txt"

tellurics_file = "tellurics_ESO.txt"

orders_file = "orders_spirou.dat"

save_file = "Simulation_HD189_ANDES.pkl"
#############################################################################

def get_rvs(t,k,p,t0):

    """
    --> Inputs:     - t:   Time vector
                    - k:   Semi-amplitude of the planet-induced RV signal on the host star
                    - p:   Planet orbital period
                    - t0:  Planet mid-transit time

    --> Outputs:    - Planet-induced RV signature for the input time values 
    """

    return - k*np.sin(2.*np.pi/p * (t-t0))

def rvp(phase,k,v):
    """
    --> Inputs:     - phase: Planet orbital phase (T-T_obs)/Porb
                    - k:     Semi-amplitude of the planet RV
                    - v:     Planet RV at mid-transit

    --> Outputs:    - Planet RV for the input orbital phases
    """
    return k*np.sin(2.*np.pi*phase)+v

def compute_transit(Rp,Rs,ip,T0,ap,Porb,ep,wp,limb_dark,uh,T_obs):

    """
    --> Inputs:     - Rp:        Planetary radius
                    - Rs:        Stellar radius (same unit as Rp)
                    - ip:        Transit inclination [deg]
                    - T0:        Mid-transit time (same unit as T_obs -- here: bjd)
                    - ap:        Semi-major-axis [Stellar radius]
                    - Porb:      Planet orbital period (same unit as T_obs)
                    - ep:        Eccentricity of the planet orbit
                    - wp:        Argument of the periapsis for the planet orbit [deg] 
                    - limb_dark: Type of limb-darkening model: "linear", "quadratic", "nonlinear" see https://lweb.cfa.harvard.edu/~lkreidberg/batman/
                    - uh:        Limb-darkening coefficients matching the type of model and in the SPIRou band (Typically H or K)
                    - T_obs:     Time vector

    --> Outputs:    - flux:      Relative transit flux (1 outside transit) 
    """
#
    params           = batman.TransitParams()
    params.rp        = Rp/Rs                       
    params.inc       = ip
    params.t0        = T0
    params.a         = ap
    params.per       = Porb
    params.ecc       = ep
    params.w         = wp         
    params.limb_dark = limb_dark
    params.u         = uh
    bat              = batman.TransitModel(params,T_obs)
    flux             = bat.light_curve(params)
    return flux


def broaden(wl,R,sigma) :

    
    c0 = 299792458.0

    #we take a kernel of size 50 000 km/s, it is exagerated but is safe

    vlim = 50000.0
    nv = 5000
    v = np.linspace(-vlim,vlim,nv)
    dv = v[1]-v[0]
    #we interpolate the model onto a regularly spaced speed array
    w0 = np.mean(wl)
    speed = c0*(w0/wl-1)
    speed_int = np.arange(0.995*np.min(speed),0.995*np.max(speed),step=dv)
    fmod = interpolate.CubicSpline(speed[::-1],R[::-1])
    mod_int = fmod(speed_int)

    

    second_conv = np.exp(-v**2/2/sigma**2) #instrumental kernel

    conv2_mod = 1/sigma/np.sqrt(2*np.pi)*2*vlim/nv*signal.oaconvolve(mod_int,second_conv,mode="same")  
    wl_int = w0/(1+speed_int/c0)


    #the model is largely oversampled and we don't want that, it is going to kill the calculation time
    diff = len(conv2_mod)/len(wl)
    spacing = max(1,round(5.*diff/6))


    return (wl_int[nv:-nv:spacing][::-1],conv2_mod[nv:-nv:spacing][::-1])