#---Notes for using Carl's codes for fitting Gaussian functions to HI emission/absorption pairs

(1) Obtain data in IDL .sav format. Need velocity arrays, emission and absorption spectrum arrays, error arrays for each.
(2) Modify cnm.pro with initial guesses for the height (hgt0), width (wid0) and center (cen0) of each Gaussian you want to fit to the absorption spectrum. Explanations of the parameters set in cnm.pro can be found in the header for gfitflex_exp.pro, written by Carl.
(3) Modify wnm.pro with initial guesses for the additional components fit to the emission profile. Explanations of the parameters set in cnm.pro can be found in the header for tbgfitflex_exp.pro, written by Carl.
(4) Modify the directory ("dir") in fit.pro where the output files will go. Set simple=1 and wfit=1 and specify the sourcename, then run fit!
