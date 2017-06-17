# Aperture-Photometry-in-DS9
A plugin python script for aperture photometry on files in DS9 with a single keypress.

## Installation
A completely default installation will clone the repository to \opt\local\bin. Customization is available and detailed below. 

The file "phot.ds9" must be placed in a directory searched by DS9 in order for single-keypress activation of the analysis script. In the default configuration, this directory is /opt/local/bin, but others are possible: see http://ds9.si.edu/doc/ref/analysis.html

Within the phot.ds9 file, you can change the installation location of the analysis script ds9_phot.py by altering the last line of the file. By default, the location is set to /opt/local/bin/ds9\_phot.py

If you choose a different installation directory, you must also alter the first line of ds9_phot.py to reflect this. The supporting file phot\_vars.py must be in the same directory as ds9\_phot.py

### Dependencies
The analysis script requires numpy, scipy, matplotlib, and photutils. The analysis script ds9\_phot.py uses Tkinter for graphical display. 

The package was built on Ubuntu 14.04 / 16.04 with package versions:
  * numpy=1.12.1
  * scipy = 0.19.0
  * matplotlib==2.0.0
  * photutils==0.3.2

After development is mostly complete, I will check compatibility with older versions on another computer.
