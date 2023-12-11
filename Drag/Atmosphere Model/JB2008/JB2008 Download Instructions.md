The JB2008 model is too large to upload to github. To compute your own versions of the density variables (i.e. for different altitude ranges/discretizations or custom space weather assumptions), please follow the following instructions.

### Installation instructions
1. Download and install SPICE Toolkit for Matlab: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
2. Set the path to the SPICE Toolkit directory in dens_prediction.m
3. Download SPICE kernels (i.e. ephemeris files) from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ and put them in the folder Data. See the links below. 

### Ephemeris files
Download the following ephemeris files and put them in the \Drag\Atmosphere Model\JB2008\JB2008_dens_field_generation\Data folder:
* de430.bsp:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
* The package comes preloaded with some additional required ephemeris files. 

### Use Instructions
1. Replace yearssolar.mat with your information and update the line "alt = linspace(200,2000,37); %km" in dens_prediction with your preferred values for minimum altitude [km]. maximum altitude [km], and number of shells.
2. Regenerate the density field using dens_prediction. This may take several hours.
3. Save the resulting .mat file to somewhere within your MATLAB path variable.