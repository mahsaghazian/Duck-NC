% CONTENTS of Elements of a Study of the Weather-within-climate variability
% and Downscaling (WEACLIM). These tools are dedicated to the analysis and
% simulation of daily time series of rainfall. Of course, several functions 
% could be used for other time series.
%
% Matlab functions written under Matlab 6. This tool needs statistic and
% signal toolboxes
%
% Version 1, April 27 2006
% suggestions, comments and bug reports to Vincent Moron (moron@cerege.fr,
% vincent@iri.columbia.edu) if you use some or all of these tools qnd if
% you need specific help.
%
% ALPHABETICAL ORDER
% anomaly - computes the anomalies relative to the seasonal cycle
% bp - barnett-Preisendorfer Canonical Correlation Analysis
% chi2 - Khi-2 test for contingeny table
% choose_area_gridded - selects a geographical window in a gridded dataset
% choose_date_daily - selects a temporal window in a dataset (daily time
% scale)
% choose_date_monthly - selects a temporal window in a dataset (monthly time scale)
% classe - ranks observations in class
% copy - copy a matrix/vector.
% corr_mc - computes a correlation and performs a random-phase test
% cv - computes cross-validation index
% decalage - shift columns of a matrix (code written by J.M. Gaillard, Dijon, France)
% distance_euclid - computes Eudlidean Distance between vectors
% dof - computes spatioo-temporal degrees of freedom (DOF)
% dof_comb - computes DOF on exhaustive combinations of columns of a matrix
% dof_comb_red - computes DOF on combinations of columns of a matrix
% ebisuzaki - create simulated random time series having the same spectral density 
% but random phases as an input time series
% filtrage - filter time series using a recursive Butterworth filter
% find_analog - search analog of a state vector in a library
% find_spell - computes statistics of dry/wet spells in time series
% fliplonlat - transform 2D matrix with first columns describing longitudes
% or latitudes into reversed format
% inverse_fisherz - performs the inverse of a fisher-Z transform
% inverse_log_odd - performs the inverse of log_odd transform
% inverse_persist - transforms back x,y into p01 and p11 (used in
% dowsncaling studies)
% knn_analog_model - DRIVER for computing knn analog analysis between observed and simulated
% fields to estimate a target observed variable. 
% lepscont -  This function computes linear error in painting space (LEPS) for
% continuous forecast (code written by N. Philippon and P. Oettli, Dijon,
% France)
% local_scaling -  This function helps to scale/calibrate a GCM time series so that the
% frequency of wet days and mean intensity of rain during wet days matches
% the observed one. 
% log_odd - performs log_odd transform
% markov_chain_1 - generates Markov chains of 0 and 1 that obeys to a Markov chain of order
% 1 (i.e. the state at t depends only on state at day at t-1).
% mixed_exp_rnd - generates simulated time series having a mixed
% exponential PDF
% mixed_expfit - computes the parameters of a mixed exponential PDF (code
% written by T. Lane, USA)
% mos_cca - DRIVER for cross validated reconstruction of a predictand field 'X' from a predictor
% field 'Y'
% nancov - computes a covariance matrix with missing values
% nandof - computes DOF with missing values
% nanstan - computes anomalies with missing values
% padd - padd matrices with values
% persist - computes x,y parameters from p01 and p11 probability transitions
% ponder - computes spatial ponderation using latitude
% prob_wet_dry - computes probability transition p00, p01, p10 and p11 in
% time series
% proj - projects from one input grid to an output grid
% replace_missing - replaces missing entries with long-term mean
% roc - computes Relative Operating Characteristics
% rocscore - computes ROC score and curves
% rps - computes Rank Probability Score
% scale_mean_var - calibrates mean and variance of a vector with the mean
% and variance of another vector
% search_mode - used in cross-validation
% select_time_area - selects temporal and spatial windows in a gridded
% dataset
% shiftgrid - shift a gridded dataset around 0° of longitude
% stan - standardizes columns of a time series
% subseas - computes seeveral quantities characterizing the interannual
% variability of seasonal and subseasonal variations of daily rainfall.
% swg_cluster - generates stochastic time series matching markov chain properties
% swg_gamma - generates stochastic time series (order 1 of markov chain)
% and gamma PDF)
% swg_mixexp - generates stochastic time series (order 1 of markov chain)
% and mixed exponential PDF)
% swg_weibull - generates stochastic time series (order 1 of markov chain)
% and weibull PDF)
% test_markov_chain - compute significance levels of probability
% transitions between types
% weather_type_model - DRIVER for computing weather types between observed and simulated
% fields to estimate a target observed variable.
%
% THEMATIC
% 
% GENERAL DRIVERS FOR ANALYSIS AND SIMULATION OF DAILY TIME SERIES
% mos_cca - DRIVER for cross validated reconstruction of a predictand field 'X' from a predictor
% field 'Y'
% local_scaling -  This function helps to scale/calibrate a GCM time series so that the
% frequency of wet days and mean intensity of rain during wet days matches
% the observed one. 
% knn_analog_model - DRIVER for computing knn analog analysis between observed and simulated
% fields to estimate a target observed variable. 
% weather_type_model - DRIVER for computing weather types between observed and simulated
% fields to estimate a target observed variable.
%
% EIGENANALYSIS - STATISTICAL ANALYSIS
% bp - barnett-Preisendorfer Canonical Correlation Analysis
% dof - computes spatioo-temporal degrees of freedom (DOF)
% dof_comb - computes DOF on exhaustive combinations of columns of a matrix
% dof_comb_red - computes DOF on combinations of columns of a matrix
% mos_cca - DRIVER for cross validated reconstruction of a predictand field 'X' from a predictor
% field 'Y'
% nancov - computes a covariance matrix with missing values
% nandof - computes DOF with missing values
%
% MANIPULATION OF GRIDDED MATRICES
% choose_area_gridded - selects a geographical window in a gridded dataset
% choose_date_daily - selects a temporal window in a dataset (daily time
% scale)
% choose_date_monthly - selects a temporal window in a dataset (monthly time scale)
% copy - copy a matrix/vector.
% decalage - shift columns of a matrix (code written by J.M. Gaillard, Dijon, France)
% fliplonlat - transform 2D matrix with first columns describing longitudes
% or latitudes into reversed format
% padd - padd matrices with values
% proj - projects from one input grid to an output grid
% scale_mean_var - calibrates mean and variance of a vector with the mean
% and variance of another vector
% search_mode - used in cross-validation
% select_time_area - selects temporal and spatial windows in a gridded
% dataset
% shiftgrid - shift a gridded dataset around 0° of longitude
%
% MISCELLANOUS
% anomaly - computes the anomalies relative to the seasonal cycle
% classe - ranks observations in class
% cv - computes cross-validation index
% ponder - computes spatial ponderation using latitude
% replace_missing - replaces missing entries with long-term mean
% scale_mean_var - calibrates mean and variance of a vector with the mean
% stan - standardizes columns of a time series
%
% STATISTICAL TESTS
% chi2 - Khi-2 test for contingeny table
% corr_mc - computes a correlation and performs a random-phase test
% ebisuzaki - create simulated random time series having the same spectral density 
% lepscont -  This function computes linear error in painting space (LEPS) for
% continuous forecast (code written by N. Philippon and P. Oettli, Dijon,
% France)
% roc - computes Relative Operating Characteristics
% rocscore - computes ROC score and curves
% rps - computes Rank Probability Score
% test_markov_chain - compute significance levels of probability
% transitions between types
%
% STOCHASTIC WEATHER GENERATORS AND ANALYSIS OF DAILY TIME SERIES
% filtrage - filter time series using a recursive Butterworth filter
% find_analog - search analog of a state vector in a library
% find_spell - computes statistics of dry/wet spells in time series
% inverse_fisherz - performs the inverse of a fisher-Z transform
% inverse_log_odd - performs the inverse of log_odd transform
% inverse_persist - transforms back x,y into p01 and p11 (used in
% dowsncaling studies)
% log_odd - performs log_odd transform
% markov_chain_1 - generates Markov chains of 0 and 1 that obeys to a Markov chain of order
% 1 (i.e. the state at t depends only on state at day at t-1).
% mixed_exp_rnd - generates simulated time series having a mixed
% exponential PDF
% mixed_expfit - computes the parameters of a mixed exponential PDF (code
% written by T. Lane, USA)
% prob_wet_dry - computes probability transition p00, p01, p10 and p11 in
% time series
% subseas - computes seeveral quantities characterizing the interannual
% variability of seasonal and subseasonal variations of daily rainfall.
% swg_cluster - generates stochastic time series matching markov chain properties
% swg_gamma - generates stochastic time series (order 1 of markov chain)
% and gamma PDF)
% swg_mixexp - generates stochastic time series (order 1 of markov chain)
% and mixed exponential PDF)
% swg_weibull - generates stochastic time series (order 1 of markov chain)
% and weibull PDF)





