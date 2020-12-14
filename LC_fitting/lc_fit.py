"""# Import curve fitting package from scipy
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv("LC_result_file.csv",sep="\t")
fit_data = pd.DataFrame(csv)
fit_data = fit_data.loc[fit_data["Flux_error"]!=0]

x_data = 0.5*(fit_data.Tmin_MJD+fit_data.Tmax_MJD)
x_error = x_data-fit_data.Tmin_MJD


plt.figure(figsize=(16,3))
plt.errorbar(x_data, fit_data.Flux, yerr=fit_data.Flux_error, xerr=x_error, fmt='.',elinewidth=None)

def fit_function(t,const_flux,flux_t0,t_rise,t_decay):
 #const_flux = 4.899271186e-8
 #flux_t0 = 1.89150529815e-07
 t0 = 56805.353
 return const_flux+(flux_t0/(np.exp((t0-t)/t_rise)+np.exp((t-t0)/t_decay)))


pars, cov = curve_fit(f=fit_function, xdata=x_data, ydata=fit_data.Flux,)
print(pars)
xfit = np.arange(56740.003,56891.003,1)
plt.plot(xfit,fit_function(xfit,*pars))

plt.savefig("fig.png")
plt.show()
"""


# Import curve fitting package from scipy
from lmfit import Model,CompositeModel
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

csv = pd.read_csv("LC_result_file.csv", sep="\t")
fit_data = pd.DataFrame(csv)
fit_data=fit_data.loc[fit_data["Flux_error"]!=0]

#print(fit_data.mean())
x_data = 0.5*(fit_data.Tmin_MJD+fit_data.Tmax_MJD)
x_error = x_data-fit_data.Tmin_MJD

def fit_function(t,t1,const_flux_1,flux_t1,t_rise_1,t_decay_1,t2,const_flux_2,flux_t2,t_rise_2,t_decay_2):
    peak1 = const_flux_1+flux_t1/(np.exp((t1-t)/t_rise_1)+np.exp((t-t1)/t_decay_1))
    peak2 = const_flux_2+flux_t2/(np.exp((t2-t)/t_rise_2)+np.exp((t-t2)/t_decay_2))
    return peak1+peak2


model = Model(fit_function,prefix="p1_")

model.set_param_hint("t1",value = 56805.303,min=-np.inf,max=np.inf,vary=True)
model.set_param_hint("t2",value = 56850.303,min=-np.inf,max=np.inf,vary=False)

model.set_param_hint("const_flux_1",value = .674549e-08,min=-np.inf,max=np.inf,vary=False)
model.set_param_hint("const_flux_2",value = 2.674549e-08,min=-np.inf,max=np.inf,vary=False)

model.set_param_hint("flux_t1",value = 2.6e-07,min=-np.inf,max=np.inf,vary=False)
model.set_param_hint("flux_t2",value = 3.26e-07,min=-np.inf,max=np.inf,vary=False)

model.set_param_hint("t_rise_1",value = 7,min=-np.inf,max=np.inf,vary=False)
model.set_param_hint("t_decay_1",value = 5,min=-np.inf,max=np.inf,vary=False)

model.set_param_hint("t_rise_2",value = 6,min=-np.inf,max=np.inf,vary=True)
model.set_param_hint("t_decay_2",value = 5,min=-np.inf,max=np.inf,vary=False)

params = model.make_params()

tx = np.linspace(56705,56950,1000)
result = model.fit(fit_data.Flux,params,t = x_data,nan_policy="propagate")

comps = result.eval_components(t = tx)

fig,ax =plt.subplots(figsize=(16,3))
plt.errorbar(x_data, fit_data.Flux, yerr=fit_data.Flux_error, xerr=x_error, fmt='.',elinewidth=None)
plt.plot(tx,comps["p1_"],color = "r")

print(result.fit_report(show_correl=False))

plt.savefig("fig.png")
