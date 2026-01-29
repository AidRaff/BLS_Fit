def spherical_hn1(n,x):
    # compute value of spherical Hankel function of the first kind of order n at x. Related to Hankel function by sqrt(pi/2x) prefactor and shift of order by + 0.5.
    import numpy as np
    import scipy as sp #necessary packages
    return np.sqrt(np.pi/(2*x))*sp.special.hankel1(n+0.5,x) #compute values(s)

def spherical_hn2(n,x):
    # compute value of spherical Hankel function of the second kind of order n at x. Related to Hankel function by sqrt(pi/2x) prefactor and shift of order by + 0.5.
    import numpy as np
    import scipy as sp #necessary packages
    return np.sqrt(np.pi/(2*x))*sp.special.hankel2(n+0.5,x) #compute values(s)

def spherical_hn1d(n,x):
    #compute value of derivative of spherical Hankel function of first kind. Recurrence relation used comes from https://www.wolframalpha.com/input?i=derivative+of+spherical+Hankel+function.
    return (-1/(2*x))*(spherical_hn1(n,x) + x*(spherical_hn1(n+1,x)-spherical_hn1(n-1,x)))

def spherical_hn2d(n,x):
    #compute value of derivative of spherical Hankel function of second kind. Recurrence relation used comes from https://www.wolframalpha.com/input?i=derivative+of+spherical+Hankel+function+of+second+kind
    return (-1/(2*x))*(spherical_hn2(n,x) + x*(spherical_hn2(n+1,x)-spherical_hn2(n-1,x)))

def RBjn(n,x):
    #compute value of Riccati-Bessel function of first kind with order n at x. Related to spherical Bessel function of first kind through multiplication by argument. 
    import scipy as sp
    return x*sp.special.spherical_jn(n,x)

def RByn(n,x):
    #compute value of Riccati-Bessel function of second kind with order n at x. Related to spherical Bessel function of second kind (a.k.a Neumann function) through multiplication by argument. 
    import scipy as sp
    return x*sp.special.spherical_yn(n,x)

def RBhn1(n,x):
    #compute value of Riccati-Bessel function of third kind with order n at x. Related to spherical Bessel function of third kind (spherical Hankel function of first kind) through multiplication by argument. 
    return x*spherical_hn1(n,x)

def RBhn2(n,x):
    #compute value of Riccati-Bessel function of fourth kind with order n at x. Related to spherical Bessel function of fourth kind (spherical Hankel function of second kind) through multiplication by argument. 
    return x*spherical_hn2(n,x)

def RBjnd(n,x):
    #compute value of derivative of Riccati-Bessel function of first kind with order n at x. 
    import scipy as sp
    return x*sp.special.spherical_jn(n,x,True) + sp.special.spherical_jn(n,x)

def RBynd(n,x):
    #compute value of derivative of Riccati-Bessel function of second kind with order n at x. 
    import scipy as sp
    return x*sp.special.spherical_yn(n,x,True) + sp.special.spherical_yn(n,x)

def RBhn1d(n,x):
    #compute value of derivative of Riccati-Bessel function of third kind with order n at x. 
    return (1/2)*(spherical_hn1(n,x) - x*(spherical_hn1(n+1,x) - spherical_hn1(n-1,x)))

def RBhn2d(n,x):
    #compute value of derivative of Riccati-Bessel function of fourth kind with order n at x. 
    return (1/2)*(spherical_hn2(n,x) - x*(spherical_hn2(n+1,x) - spherical_hn2(n-1,x)))

def Mie_an_old(n, m, x):
    #compute Mie a coefficient of order n for a particle of relative refractive index m and size parameter x.
    return ((m*RBjn(n,m*x)*RBjnd(n,x))-(RBjn(n,x)*RBjnd(n,m*x)))/((m*RBjn(n,m*x)*RBhn1d(n,x))-(RBhn1(n,x)*RBjnd(n,m*x)))

def Mie_bn_old(n, m, x):
    #compute Mie b coefficient of order n for a particle of relative refractive index m and size parameter x.
    return ((RBjn(n,m*x)*RBjnd(n,x))-(m*RBjn(n,x)*RBjnd(n,m*x)))/((RBjn(n,m*x)*RBhn1d(n,x))-(m*RBhn1(n,x)*RBjnd(n,m*x)))

def RB_ratio(n,mx):
    import numpy as np
    r = 1/np.tan(mx)
    for i0 in range(n):
        r = 1/(((2*i0+1)/mx) - r)
    return r

def Mie_an(n, m, x):
    r = RB_ratio(n,m*x)
    return (((r/m) + ((n*(1-(1/m**2)))/x))*RBjn(n,x) - RBjn(n-1,x))/(((r/m) + ((n*(1-(1/m**2)))/x))*RBhn1(n,x) - RBhn1(n-1,x))

def Mie_bn(n, m, x):
    r = RB_ratio(n,m*x)
    return (r*m*RBjn(n,x) - RBjn(n-1,x))/(r*m*RBhn1(n,x) - RBhn1(n-1,x))

def Mie_pin(n, theta):
    #angular function pi of order 1-n evaluated at theta. Defined and calculated as per Bohren and Huffman P 94-95. n must be >= 1.
    import numpy as np
    if n == 1:
        return np.array([[np.ones_like(theta)]]) #pi_1 = 1 by definition
    else: 
        pin = np.zeros((n+1, len(theta))) #space for pi_n values
        pin[0, :] = np.zeros_like(theta) #row for pi_0
        pin[1, :] = np.ones_like(theta) #row for pi_1
        for i1 in range(2, n+1):
            pin[i1,:] = (((2*i1)-1)/(i1-1))*np.cos(theta)*pin[i1-1,:] - (i1/(i1-1))*pin[i1-2,:] #calculate higher order from recurrence relation
    return pin[1:,:] #order 0 not relevant for Mie theory so don't output

def Mie_taun(n, theta):
    #angular function tau of order 1-n evaluated at theta. Defined and calculated as per Bohren and Huffman P 94-95. n must be >= 1.
    import numpy as np
    pin = Mie_pin(n, theta) #first calculate pi_n for each n at each theta
    taun = np.zeros_like(pin) #space for values of tau_n
    for i1 in range(n):
        if i1 == 0:
            taun[i1,:] = np.cos(theta)*pin[i1,:]
        else:
            taun[i1,:] = (i1+1)*np.cos(theta)*pin[i1,:] - (i1+2)*pin[i1-1,:] #compute tau_n for each value of n (indices are a little confusing)
    return taun

def n_osc(nu, A, nu0):
    #Calculate real part of refractive index at wavenumber nu with oscillator parameters amplitude A and central frequency nu0. Uses approximate experession from https://doi.org/10.1029/2019gl084568
    import numpy as np
    return 1 + (2/np.pi)*((A*nu0)/(nu0**2 - nu**2))

def Pearson_CC(vec1, vec2):
    #calculates the Pearson correlation coefficient (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between two vectors. Lengths must match. 
    import numpy as np
    N = len(vec1) #number of samples
    sum1 = np.einsum('i->',vec1)
    sum2 = np.einsum('i->',vec2)
    return ((N*np.dot(vec1,vec2)) - (sum1*sum2))/((np.sqrt((N*np.dot(vec1,vec1))-((sum1)**2)))*(np.sqrt((N*np.dot(vec2,vec2))-((sum2)**2))))

def parabolic_interpolation(x, y):
    # Performs parabolic interpolation on the points x and y. Each array must contain three entries, corresponding to the x and y coordinates of points surrounding a maximum/minimum. Output is the x cooridnate of the extremum.
    d1 = y[1] - y[0]
    d2 = y[2] - 2*y[1] + y[0]
    return x[2] - ((x[2]-x[1])*((d1/d2)+1.5))

def parabolic_interpolation_xy(x, y):
    # Performs parabolic interpolation on the points x and y. Each array must contain three entries, corresponding to the x and y coordinates of points surrounding a maximum/minimum. Output is the x cooridnate of the extremum.
    import numpy as np
    px = (y[1]-y[0])/(y[2]-2*y[1]+y[0])
    py = (y[0]-y[2])/(y[2]-2*y[1]+y[0])
    xmax = x[2] - ((x[2]-x[1])*(px+1.5))
    ymax = y[1] - 0.25*(y[0]-y[2])*(py/2)
    return np.array([xmax,ymax])

def colour_subset(n, cm):
    #generates array of nx4 evenly spaced rgb and transparency values from the colormap cm
    import matplotlib as mpl
    subset = mpl.colormaps[cm].resampled(n) #resample colormap at evenly spaced values
    return subset(range(n)) #return as array

def vector_normalise(vec):
    #normalises a vector to between zero and one
    zerod = vec - vec.min() #subtract minimum value from each entry so minimum of vector is zero
    return zerod/zerod.max() #divide zerod vector by maximum to get maximum value of one

def Pearson_CC(vec1, vec2):
    #calculates the Pearson correlation coefficient (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between two vectors. Lengths must match. 
    import numpy as np
    N = len(vec1) #number of samples
    sum1 = np.einsum('i->',vec1)
    sum2 = np.einsum('i->',vec2)
    return ((N*np.dot(vec1,vec2)) - (sum1*sum2))/((np.sqrt((N*np.dot(vec1,vec1))-((sum1)**2)))*(np.sqrt((N*np.dot(vec2,vec2))-((sum2)**2))))

def generate_spectrum(r, wl, m, theta_i, phi_i, NA):
    # Generate scattering spectra for a particle of radius r across the wavelength range wl with refractive index n + ik. Currently only accepts multiple entries for radius and wavelength (not multiple RIs, will accept single entry for radius). 
    # Light is incident at (theta_i, phi_i) in spherical coordinates (see https://doi.org/10.1016/j.jqsrt.2025.1097081 for geometry) and collected with an objective of numerical aperture NA.
    import numpy as np
    if np.size(m) == 1:
        m = m*np.ones_like(wl) #bug fix for if only one value of m input
    wn = (2*np.pi)/wl #wavenumber of incident light
    x = (2*np.pi*r)/wl #size parameter of each wavelength
    omax = np.round(x + (4*(x**(1/3))) + 2).astype(int) #maximum order to which sums need to be evaluated for each size parameter
    o = np.linspace(1, omax.max(), omax.max()) #list of orders used to calculate prefactors, coefficients and angular functions
    npts = int(np.round((0.75*NA*omax.max())) + 18) #number of points to use for numerical integration
    if np.mod(npts,2) == 1:
        npts += 1
    pf = ((2*o)+1)/(o*(o+1)) #prefactors for eventual sum
    vi = np.transpose(np.array([np.sin(theta_i)*np.cos(phi_i), np.sin(theta_i)*np.sin(phi_i), np.cos(theta_i)])) #incident wavevector
    an = np.zeros((omax.max(), len(x)), dtype=complex)
    bn = np.zeros_like(an)
    for i1 in range(omax.max()):
        docalc = o[i1] <= omax
        an[i1,docalc] = pf[i1]*Mie_an(int(o[i1]), m[docalc], x[docalc])
        bn[i1,docalc] = pf[i1]*Mie_bn(int(o[i1]), m[docalc], x[docalc]) #calculate given order of coefficients
    theta_s = np.linspace(0, np.arcsin(NA), npts)#range of scattered theta to look at (angle between scattered wavevector and detector axis)
    phi_s = np.linspace(0, 2*np.pi, npts) #range of scattered phi to use (angle between scattered wavevector and x-axis)
    theta_c = np.zeros((len(theta_s), len(phi_s)))
    for i1 in range(len(theta_s)):
        for i2 in range(len(phi_s)):
            vs = np.transpose(np.array([np.sin(theta_s[i1])*np.cos(phi_s[i2]), np.sin(theta_s[i1])*np.sin(phi_s[i2]), np.cos(theta_s[i1])])) #scattering vector components (minus k, which is not needed for collection angle calculation)
            theta_c[i1,i2] = np.arccos(np.dot(vi,vs)) #collection angle for scattered wavevector
    pin = np.zeros((len(theta_s), len(phi_s), omax.max()))
    taun = np.zeros_like(pin) #space for values of pi_n and tau_n (independent of size parameter except for summation maximum, so can be evaluated before loop and called on later)
    for i1 in range(len(theta_s)):
        pn = Mie_pin(omax.max(), theta_c[i1,:])
        tn = Mie_taun(omax.max(), theta_c[i1,:]) #calculate necessary range of pi_n and tau_n across collection angle range
        pin[i1,:,:] = np.transpose(pn)
        taun[i1,:,:] = np.transpose(tn) #put into matrices
    spectrum = np.zeros_like(wl) #space for spectrum at given radius
    for i1 in range(len(wl)):
        if np.size(pin[:,:,:omax[i1]]) < 65536:
            S1 = np.dot(pin[:,:,:omax[i1]],an[:omax[i1],i1]) + np.dot(taun[:,:,:omax[i1]],bn[:omax[i1],i1])
            S2 = np.dot(taun[:,:,:omax[i1]],an[:omax[i1],i1]) + np.dot(pin[:,:,:omax[i1]],bn[:omax[i1],i1])
        else:
            S1 = complex(0,0)
            S2 = complex(0,0) #initial values for summation
            for i2 in range(omax[i1]):
                S1 += (an[i2,i1]*pin[:,:,i2]) + (bn[i2,i1]*taun[:,:,i2])
                S2 += (an[i2,i1]*taun[:,:,i2]) + (bn[i2,i1]*pin[:,:,i2]) #calculate summation terms
        S11 = (1/(2*(wn[i1]**2)))*((np.abs(S2)**2) + (np.abs(S1)**2)) #calculate S11 as function of theta and phi
        theta_int_terms = np.zeros_like(S11)
        for i2 in range(len(phi_s)):
            theta_int_terms[:,i2] = np.sin(theta_s)*S11[:,i2] #calculate terms for theta integration (multiply S11 terms by sin(theta))
        phi_int_terms = ((theta_s[1] - theta_s[0])/2)*(theta_int_terms[0,:] + theta_int_terms[-1,:] + 2*np.einsum('ij->j',theta_int_terms[1:-1,:])) #perform theta integration using Simpson's 1/3 rule for regular spaced points
        spectrum[i1] = ((phi_s[1] - phi_s[0])/2)*(phi_int_terms[0] + phi_int_terms[-1] + 2*np.einsum('i->', phi_int_terms[1:-1])) #perform phi integration using Simpson's 1/3 rule and put into vector
    return spectrum

def generate_spectra_x(x, m, theta_i, phi_i, NA):
    # Generates spectra as above but with size parameter as an input.
    import numpy as np
    omax = np.round(x + (4*(x**(1/3))) + 2).astype(int) #maximum order to which sums need to be evaluated for each size parameter
    o = np.linspace(1, np.max(omax), np.max(omax)) #list of orders used to calculate prefactors, coefficients and angular functions
    npts = int(np.round(0.75*NA*omax.max()) + 18) #number of pts to use for numerical integration
    if np.mod(npts,2) == 1:
        npts += 1
    pf = ((2*o)+1)/(o*(o+1)) #prefactors for eventual sum
    vi = np.transpose(np.array([np.sin(theta_i)*np.cos(phi_i), np.sin(theta_i)*np.sin(phi_i), np.cos(theta_i)])) #incident wavevector
    an = np.zeros((omax.max(), len(x)), dtype=complex)
    bn = np.zeros_like(an)
    for i1 in range(omax.max()):
        docalc = o[i1] <= omax
        an[i1,docalc] = pf[i1]*Mie_an(int(o[i1]), m, x[docalc])
        bn[i1,docalc] = pf[i1]*Mie_bn(int(o[i1]), m, x[docalc]) #calculate given order of coefficients
    theta_s = np.linspace(0, np.arcsin(NA), npts)#range of scattered theta to look at (angle between scattered wavevector and detector axis)
    phi_s = np.linspace(0, 2*np.pi, npts) #range of scattered phi to use (angle between scattered wavevector and x-axis)
    theta_c = np.zeros((len(theta_s), len(phi_s)))
    for i1 in range(len(theta_s)):
        for i2 in range(len(phi_s)):
            vs = np.transpose(np.array([np.sin(theta_s[i1])*np.cos(phi_s[i2]), np.sin(theta_s[i1])*np.sin(phi_s[i2]), np.cos(theta_s[i1])])) #scattering vector components (minus k, which is not needed for collection angle calculation)
            theta_c[i1,i2] = np.arccos(np.dot(vi,vs)) #collection angle for scattered wavevector
    pin = np.zeros((len(theta_s), len(phi_s), omax.max()))
    taun = np.zeros_like(pin) #space for values of pi_n and tau_n (independent of size parameter except for summation maximum, so can be evaluated before loop and called on later)
    for i1 in range(len(theta_s)):
        pn = Mie_pin(omax.max(), theta_c[i1,:])
        tn = Mie_taun(omax.max(), theta_c[i1,:]) #calculate necessary range of pi_n and tau_n across collection angle range
        pin[i1,:,:] = np.transpose(pn)
        taun[i1,:,:] = np.transpose(tn) #put into matrices
    spectrum = np.zeros_like(x) #space for spectrum at given radius
    for i1 in range(len(x)):
        if np.size(pin[:,:,:omax[i1]]) < 65536:
            S1 = np.dot(pin[:,:,:omax[i1]],an[:omax[i1],i1]) + np.dot(taun[:,:,:omax[i1]],bn[:omax[i1],i1])
            S2 = np.dot(taun[:,:,:omax[i1]],an[:omax[i1],i1]) + np.dot(pin[:,:,:omax[i1]],bn[:omax[i1],i1])
        else:
            S1 = complex(0,0)
            S2 = complex(0,0) #initial values for summation
            for i2 in range(omax[i1]):
                S1 += (an[i2,i1]*pin[:,:,i2]) + (bn[i2,i1]*taun[:,:,i2])
                S2 += (an[i2,i1]*taun[:,:,i2]) + (bn[i2,i1]*pin[:,:,i2]) #calculate summation terms
        S11 = (1/2)*((np.abs(S2)**2) + (np.abs(S1)**2)) #calculate S11 as function of theta and phi
        theta_int_terms = np.zeros_like(S11)
        for i2 in range(len(phi_s)):
            theta_int_terms[:,i2] = np.sin(theta_s)*S11[:,i2] #calculate terms for theta integration (multiply S11 terms by sin(theta))
        phi_int_terms = ((theta_s[1] - theta_s[0])/2)*(theta_int_terms[0,:] + theta_int_terms[-1,:] + 2*np.einsum('ij->j',theta_int_terms[1:-1,:])) #perform theta integration using Simpson's 1/3 rule for regular spaced points
        spectrum[i1] = ((phi_s[1] - phi_s[0])/2)*(phi_int_terms[0] + phi_int_terms[-1] + 2*np.einsum('i->', phi_int_terms[1:-1])) #perform phi integration using Simpson's 1/3 rule and put into vector
    return spectrum

def approximate_fit(spec_exp, spec_the, wranges, wpts, rpts, wl, x, n, ig, lb, ub):
    #find approximate radius and oscillator parameters for the experimental spectrum spec_exp by comparing to theoretical spectra spec_the. Correlations between the two are calculated over the wavelength ranges
    # wranges, with corresponding midpoints wpts (used for fitting). wl is the wavelength points of the experimental spectra (theoretical spectra are interpolated to these in size parameter space using each radius
    # in rpts). ig, lb and ub are the initial guess, lower bound, and upper bound for the effective oscillator fit. x is the size parameter grid on which the theoretical spectra are calculated, and n the RI of each.
    import numpy as np
    import scipy as sp
    corr = np.zeros((len(spec_the), len(rpts), len(wranges))) #space for correlations
    for i1, sc in enumerate(spec_the):
        for i2, rg in enumerate(rpts):
            xint = ((wl/(2*np.pi))**2)*np.interp(2*np.pi*rg/wl, x, sc) #interpolate calculated spectra to x values produced by i2th value of radius
            corr[i1,i2,:] = [Pearson_CC(spec_exp[wr], xint[wr]) for wr in wranges] #calculate correlation coefficient for spectrum in each wavelength range
    cmaxima = corr.max(axis=0)
    nmaxima = np.zeros_like(cmaxima, dtype=int)
    for i1 in range(len(rpts)):
        nmaxima[i1,:] = np.squeeze([np.where(corr[:,i1,i2] == corr[:,i1,i2].max())[0] for i2 in range(len(wranges))])
    valid = nmaxima[:,0] == nmaxima.max(axis=1)
    corrtot = np.einsum('ij->i', cmaxima) #sum maximum correlations in each range for each radius
    maxtot = corrtot == corrtot[valid].max()
    nmaxind = np.squeeze(nmaxima[corrtot == corrtot[valid].max()])
    wpts = wpts[(nmaxind != 0) & (nmaxind != len(n)-1)]
    nmaxind = nmaxind[(nmaxind != 0) & (nmaxind != len(n)-1)] #get rid of any maxima at extrema of RI range (causes next line to fail)
    if len(nmaxind) > 0:
        npts = np.squeeze([parabolic_interpolation(n[nmi-1:nmi+2], corr[nmi-1:nmi+2,maxtot,i1]) for i1, nmi in enumerate(nmaxind)])
        res, _ = sp.optimize.curve_fit(n_osc, 1/wpts, npts, p0=ig, bounds=(lb, ub), max_nfev=1e5) #fit constant RI points to oscillator model to get rough oscillator parameters
    else:
        res = [0,0] #if no points (occasionally happens) output 0 for oscillator parameters
    return np.array([np.squeeze(rpts[maxtot]), res[0], res[1]])