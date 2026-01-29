BLS_Fit README

SUMMARY

BLS_Fit is a Python code used to fit broadband light scattering (BLS) spectra of spherical homogeneous particles to determine the particle radius and wavelength-dependent refractive index. Full details of the retrieval method can be found in the associated publication (link added upon publication). A brief summary is given below.

Users input several parameters at the start of the code:

rmin, rmax, rinc: minimum (rmin) and maximum (rmax) values of radius to include in approximate parameter search, and increment (rinc) to span this space in (10 nm recommended, but feel free to experiment). Let intuition guide the choice of rmin and rmax.

nmin, nmax, ninc: same but for refractive index values. Again, let intuition guide the choice of limits (e.g. typical values for aqueous inorganic particles would be 1.35-1.5), recommended increment is 0.005, but again feel free to experiment. 

show_examples: option to generate extra graphs at the end of the approximate and precise retrieval steps which show the application of these methods to a single measured spectrum. Intended to offer a graphical illustration of the method. Set True to show, False to turn off.

theta_i, phi_i, NA: theta_i and phi_i specify the direction of the light incident on the particle and NA the angular range over which scattered light is collected. See the associated publication for a diagram of the experimental geometry. 

spec_range_lims: reduced spectral ranges to use in approximate fitting step. Should be list of len(2) numpy arrays, with the entries in each being lower and upper limits of the windows, respectively. Discussion and guidance on how to choose these is found in the associated publication. However, you need a minimum of three windows and all should lie within the spectral range of your measurements. 

ncores: number of cores to use in parallel computations. There are two points which use parallel computing to increase computational speed. If you have only a single core available, simply set this to 1 and everything will still run fine. 

The examples included here all include some code that generates the spectra to act as data. Obviously, when you try to apply this to your measurements the code generating this "data" should be replaced with code for importing your own data. At the end of your import procedure, you need three things for the code to run:
1) A numpy array called wl which contains the wavelengths at which spectral intensities are measured.
2) A list called spectra, where each list entry is an array containing the intensities for a given frame. Each list entry should be the same length as wl.
3) A numpy array called frameno, which should be the same length as spectra and simply numbers the frames for plotting purposes. 

With the above parameters set and the necessary data imported the code should then run without further intervention. Full details of the process are described in the associated publication, but the process essentially breaks the analysis into two stages. The first step determines rough particle parameters by generating spectra with a constant refractive index for every combination of r and n specified by the user at the start. Each measured spectrum is broken into the discrete wavelength ranges specified in windows at the start of the notebook, and correlations between the measured and constant-n spectra in each window are used to generate approximate guesses for radius and the two parameters required to describe the wavelength-dependent refractive index, A and nu_0. These approximate parameters are then used to generate synthetic spectra with wavelength-dependent refractive indices for use in the second stage of the process. The second stage calculates correlations between measured spectra and the wavelength-dependent RI spectra and interpolates between those with the highest correlations to determine the precise fitting parameters r_fit, A_fit and nu_fit. These are then output as columns in a .txt file.

REQUIREMENTS

In order to run this code, you must have Python 3 installed on your computer. Those new to Python are recommended to install Anaconda (https://www.anaconda.com/download) as this is, in the author's experience, the easiest way to get started with Python. Code was written in Python 3.11.5, and has been tested with Python 3.13.9, but should run with any version of Python 3.

In the main folder, in addition to this README.txt file, you will find a license file (LICENSE.txt), a Python file containing functions necessary for running the code (BLS_Fit_Functions.py), and a minimal working example in the form of a Jupyter notebook, BLS_Fit.ipynb. The minimal working example generates some synthetic data with a constant radius and refractive index and demonstrates how the fitting procedure is applied to this simple case. The notebook format is used to integrate diagnostic graphs throughout the process so users can see how the analysis proceeds and identify any potential problems that need their attention. The folder Examples will contain more complex examples, including the two cases shown in the associated publication: evaporation from a single component particle, and the RH response of a hygroscopic particle. These demonstrate how the method can be applied to more complex and realistic scenarios. 

In addition to the BLS_Fit_Functions.py file, which should be kept in the same folder as the notebook for purposes of running the code, users are also required to have several other packages installed. All should be part of any standard Python IDE, but can be downloaded individually too. The required packages are:

1) NumPy (A package providing several basic computational functions, downloadable here: https://numpy.org/)
2) SciPy (A package providing other functions specifically for scientific computing, downloadable here: https://scipy.org/)
3) Matplotlib (A package used for plotting graphs, downloadable here: https://matplotlib.org/)
4) Time (A package for timing code execution, part of the Python Standard Library. More info here: https://docs.python.org/3/library/time.html)
5) Joblib (A package used for running calculations in parallel, downloadable here: https://joblib.readthedocs.io/en/stable/)

AUTHORS

This code is written and maintained by Aidan Rafferty. He can be contacted at aidan.rafferty@mail.mcgill.ca. 

COPYRIGHT AND LICENSING

This work is licensed under the CC BY 4.0 License (more information in LICENSE.txt or at https://creativecommons.org/licenses/by/4.0/).