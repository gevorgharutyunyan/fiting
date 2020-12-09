# Import curve fitting package from scipy

from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv("LC_result_file.csv",sep="\t")
fit_data = pd.DataFrame(csv)

x_data = 0.5*(fit_data.Tmin_MJD+fit_data.Tmax_MJD)
x_error = x_data-fit_data.Tmin_MJD

plt.figure(figsize=(16,3))
plt.errorbar(x_data, fit_data.Flux, yerr=fit_data.Flux_error, xerr=x_error, fmt='.',elinewidth=None)

def fit_function(t,t_rise,t_decay):
    const_flux = 5.037271186e-8
    flux_t0 = 1.89150529815e-07
    t0 = 56805.353
    return const_flux+(flux_t0/(np.exp((t0-t)/t_rise)+np.exp((t-t0)/t_decay)))


pars, cov = curve_fit(f=fit_function, xdata=x_data, ydata=fit_data.Flux)
print(pars)
xfit = np.arange(56740.003,56891.003,1)
plt.plot(xfit,fit_function(xfit,*pars))


plt.show()
"""# Function to calculate the Gaussian with constants a, b, and c

# Generate dummy dataset
x_dummy = np.linspace(start=-10, stop=10, num=100)
y_dummy = fit_function(x_dummy, 8, -1, 3)
# Add noise from a Gaussian distribution
noise = 0.5*np.random.normal(size=y_dummy.size)
y_dummy = y_dummy + noise


# Fit the dummy Gaussian data

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevs = np.sqrt(np.diag(cov))
# Calculate the residuals
res = y_dummy - power_law(x_dummy, *pars)"""