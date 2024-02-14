# Sequence generator


## Activities

The activity generation simulates an experimental data form some real assay with noise and bias.

For each sequence the algoritm is folowing
* Generate "ideal" activity - standard Gauss-distributed value centered at 0 and with dispersion 1.
* Generate "noise" - Gauss-distributed value with a randomly chosen center in the interval [-2,2] 
  and randomly chose dispersion in the interval [0,2]
* Calculate real_actilivy = (ideal_activity + noise*noise_level)/(1+noise_level)
  Noize_level is the ration of signal/noize.

* Add scaling and systemstic error of the assay assay_scale is some ran
  * assay_scale - scaling parameter - default is a random value in [0-10].
  * assay_average - the mean of assay default is a random value in [0-10].
* measured_activity = real_activity * assay_scale + assay_avarage.
