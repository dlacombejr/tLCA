Auto Gaussian & Gabor fits
---
Functions to fit a 1D Gaussian to a curve and a 2D Gaussian or Gabor to a surface. The routines are automatic in the sense that they do not require the specification of starting guesses for the model parameters. This is done by evaluating the quality of fit for many different choices of parameters then refining the most promising set of params through least-squares (exhaustive search followed by refinement). 

All functions support 2 methods for computing error bars on the parameters: bootstrapping and MCMC.

autoGaussianSurf(xi,yi,zi) fits a 2D Gaussian to a surface, defined as:

zi = a*exp(-((xi-x0).^2/2/sigmax^2 + (yi-y0).^2/2/sigmay^2)) + b

autoGaborSurf(xi,yi,zi) fits a Gabor, defined as:

zi = a*exp(-(xip,.^2+yip.^2)/2/sigma^2)*cos(2*pi*xip/lambda + phase) + b

Where:
xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta);
yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);

The Gabor fit calls autoGaussianSurf internally, using the fact that the absolute value of a Gabor in the Fourier domain is a Gaussian.

autoGaussianCurve(xi,zi) fits a 1D Gaussian to a curve.

Author: Patrick Mineault
patrick DOT mineault AT gmail DOT com

History:
21-10-2011 - Added 1D Gaussian curve fitting
           - Added MCMC and bootstrapping-based error bars
17-09-2011 - Added MCMC method for error bars (Metropolis-Hastings)
01-08-2011 - Removed Gibbs sampling version of Gaussian surface fit (was unreliable).
           - Added Gabor fitting
           - Changed function names
           - Uses better limits for sigmax, sigmay for Gaussian fit
08-06-2011 - Included C file for Gibbs sampling version
           - Added some basic and not entirely reliable convergence checks for Gibbs sampling
18-05-2011 - Initial release
