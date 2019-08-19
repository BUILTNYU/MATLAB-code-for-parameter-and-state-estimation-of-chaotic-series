# MATLAB-code-for-parameter-and-state-estimation-of-chaotic-series
This is code that I wrote in 2011 for jointly estimating the parameters and states of a chaotic series. The algorithm is from: 

Nakamura, T., Hirata, Y., Judd, K., Kilminster, D., 2007. Improved parameter estimation from noisy time series for nonlinear dynamical systems. International Journal of Bifurcation and Chaos 17(5), 1741-1752. 

I modified their method to have a smoother descent using Method of Successive Averages. The purpose was to apply this method to estimate chaotic arrivals in a chaos-driven queue, which is both deterministic and unpredictable, allowing one to analyze transient properties and propagate complex queue patterns over a network. This work is in: 

Chow, J.Y.J., 2013. On observable chaotic maps for queueing analysis. Transportation Research Record, accepted for publication.
