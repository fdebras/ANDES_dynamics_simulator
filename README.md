# ANDES_dynamics_simulator
A simple simulator for ANDES data

If you don't want to get bothered and just load the data, I provide a Simulation_HD189_ANDES.pkl file which has a SNR 10 times that of SPIRou (hence same exposure time at first order) which you can read by : 

with open("Simulation_HD189_ANDES.pkl",'rb') as specfile:
    orders,wl,data,blaze,tellurics,T_obs,phase,window,berv,vstar,airmass,SN = pickle.load(specfile)

It is accessible here: https://drive.google.com/file/d/1mAQmr4YbZQ1u6rs5T6AN9KPMlV9XsT3D/view?usp=sharing 

Using the code: 

You will need the SPIRou data file, that you can download here : https://drive.google.com/file/d/11sOBNFa0Yf6FyoztBmjQfXTq8YvI4OgO/view?usp=sharing

The code is called "create_ANDES_data.py" and relies on the parameters and fuctions defined in "parameters_and_functions.py". What I decided to do is to scale everything from a SPIRou transit of HD 189733b obtained in June 2020 as it was the easiest in the first place, and I am now reading the docs we got from Nicoletta to update that more specifically to what we expect of ANDES. I have not considered the visible for now because I think that most of the dynamics we are interested in is in the infrared but I will also need to do it with the ANDES documents. 

After reading the data and interpolating the star and planet spectra (the star data are its flux and the planet data are the planet radius in cm here), I create the obersvation time ( line 45) based on the requested SNR: as the ELT is 10 times bigger, we have 10 times more SNR in a same exposition time (100* flux). Hence there is a scaling parameter: if you put it to 10, then there will be the same number of exposures than a SPIRou transit but with 10 times more SNR, if you put it to 1 you will have 100 times more exposures than a SPIRou transit with the same SNR. From this array of observation time we can then define the phase, window, airmass and stellar/planetary RV during transit. 

Starting from line 74, we then cut the data into chunk of wavelength following the SPIRou orders. We do this at at resolution 5 times that of ANDES. We create the noise (assuming white noise) following the SNR distribution of SPIRou, as well as a blaze following the blaze function of SPIRou order per order. 

Line 117 we then interpolate the stellar spectra at the appropriate speed. I then add the possibility to "consider flux variations": with SPIRou, we cannot know how much flux of the star entered the fiber because the positioning of the fiber varies with time and the monitoring is not precise enough. Hence I added the possibility that each spectra has a different mean flux, as in SPIRou, and inspired from the variance observed in the SPIRou sequence of june 2020.

We then add the planet line 129. We should add tellurics afterwards but I mentioned my issue already. 

And we can now broaden the model at the instrumental resolution, before averaging this broadening over a pixel. The avergaing here is very naive. Note that I made the same assumption than in SPIRou : a pixel is half the resolution size. 

Then this is it, we just have to save everything up. 
