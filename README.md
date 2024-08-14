# Early indicator for bounded noise model
Numerical approximation of derivative of extremal maps for one-dimensional random difference equation with bounded noise. This provides an early indicator for critical transitions.

File 'lambda_estimate.m'

The first part of the code simulate $n$ number of data with different number of noise realisations on the variable 'sample_size'. 

In the second part, we approximate the derivative of the extremal map by constructing a normalised histogram data with number of bins 'BinNo' and portion of data 'data_portion'. Our first quadratic method in the paper, we set 'brutequad = 1', while for the improved method, we set in addtion 'brutequad11 = 1'. 

The third part of the code uses the interval method from the literature - C. Kuehn, G. Malavolta, and M. Rasmussen, Early-warning signals for bifurcations in random dynamical systems with bounded noise, Journal of Mathematical Analysis and Applications, 464 (2018), pp. 58-77.

The final part is the plotting of graphs.

File 'saddle_discrete.m' discribe different maps used for simulations, including the hypertangent map or linear map with different types of noise.
