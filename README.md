# growth-rates
Code for inferring time-varying epidemic growth rates (r_t). Two methods are used:

i) The Wallinga-Lipsitch method which converts the time-varying reproduction number (R_t) estimated from EpiFilter (see https://github.com/kpzoo/EpiFilter) into growth rates under a gamma distributed serial interval assumption.
ii) A generic smoothing approach defining the growth rate as the smoothed derivative of the logarithm of the incidence curve - based on the theory of Savitsky Golay filtering, which includes many approaches such as least squares spline fitting.

These analyses reproduce the Figures from: Parag KV, Thompson RN, Donnelly CA. Are epidemic growth rates more informative than reproduction numbers? medRxiv. 2021; doi:10.1101/2021.04.15.21255565. See https://rss.org.uk/RSS/media/File-library/Publications/Special%20topic%20meeting/Final_Manuscript_Parag_Thompson_Donnelly.pdf
for the paper, which will appear as part of a JRSSA special issue. A talk on this work is https://www.youtube.com/watch?v=GTKu4fjA0Fc&ab_channel=RoyalStatSoc.

The Matlab code consists of growthSim.m, which computes the comparisons of growth rate estimates under a known serial interval distribution, and growthSimErr.m, which examines how estimates of R_t and r_t are differently sensitive to misspecification of that distribution.
