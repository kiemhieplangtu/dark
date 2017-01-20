*** How to fit emission/absorption spectra *** 
3C286, one of the 21-SPONGE sources (Figure 3 of Murray+15), is used here as an example. 
(1) Fit the absorption spectrum first by providing initial guesses in "cnm.pro".
--- This can be done by eye-examining the absorption spectrum, e.g., 
IDL> restore, /ver, '3C286.sav'
IDL> plot, vlsr[where((vlsr gt -60.0) and (vlsr lt 30.0))], spec1[where((vlsr gt -60.0) and (vlsr lt 30.0))], yrange=[0.99,1.002], ystyle=1
--- Goodness-of-fit can be examined by comparing "rms" and "sigma", e.g., 
IDL> restore, /ver, '3C286.sav'
IDL> print, rmsprint
     0.000696729
IDL> restore, /ver, '3C286_ABS_params.sav'
IDL> print, sigma
     0.000806256
    CNM components can be added until "sigma" becomes smaller than or comparable to "rms". 
(2) Then fit the emission spectrum by adding WNM components in "wnm.pro". 
--- Like the absorption fitting case, initial guesses can be provided by eye-examining the emission spectrum, e.g., 
IDL> plot, vlsrem[where((vlsrem gt -60.0) and (vlsrem lt 30.0))], em_spec[where((vlsrem gt -60.0) and (vlsrem lt 30.0))]
--- Like the absorption fitting case, goodness-of-fit can be examined by 
IDL> restore, /ver, '3C286.sav'
IDL> emsigma_sub1 = em_spec[where((vlsrem gt -80.0) and (vlsrem lt -50.0))]
IDL> emsigma_sub2 = em_spec[where((vlsrem gt 50.0) and (vlsrem lt 80))]
IDL> emsigma_sub = [emsigma_sub1, emsigma_sub2]
IDL> print, stddev(emsigma_sub)
     0.0710359
IDL> restore, /ver, '3C286_EM_params.sav'
IDL> print, sigmaw
     0.18421481     
    WNM components can be added until "sigmaw" becomes reasonably close to the 1-sigma uncertainty of the emission spectrum. 

Note: 
(1) "vlsr" and "vlsrem" arrays have different sizes because they are from different telescopes (VLA and Arecibo). 
(2) Currently, "fit.pro" does not take into account different spectral resolutions for Gaussian decomposition. 
(3) Claire's advice on Gaussian decomposition is: try to add minimum WNM components!
Otherwise, the derived spin temperature can be unrealistically low. 
Remember Tb = Ts(1-exp(-tau)) and tau is provided from the absorption fitting.

Last updated: 09/04/2015  
