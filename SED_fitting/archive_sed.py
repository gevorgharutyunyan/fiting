import matplotlib.pyplot as plt
import pandas as pd
from Constants import*
from scipy.integrate import quad

saqo_txt = pd.read_csv("saqo_check.txt", header=None, delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False)
saqo = pd.DataFrame(saqo_txt)



archiv_data = pd.read_csv("PKS1430178.txt", header=None, delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False)
archiv_sed = pd.DataFrame(archiv_data)

flare_data = pd.read_csv("flare_pks_1430_178_sed_data.csv", delim_whitespace=True,
                         error_bad_lines=False, warn_bad_lines=False)
flare = pd.DataFrame(flare_data)
flare=flare[flare.Flux_Error !=0]

quit_data = pd.read_csv("pks_1430_178_sed_data.csv", delim_whitespace=True,
                        error_bad_lines=False, warn_bad_lines=False)
quit = pd.DataFrame(quit_data)
quit=quit[quit.Flux_Error !=0]

fl_err1 = flare.Energy_Center-flare.Energy_Min
fl_err2 = flare.Energy_Max-flare.Energy_Center

qt_err1 = quit.Energy_Center-quit.Energy_Min
qt_err2 = quit.Energy_Max-quit.Energy_Center


fig, ax = plt.subplots(figsize=(14, 5))

#for check
ax.plot(saqo[0],saqo[1],color='y',)
# Archived data from asi.it
#ax.errorbar(hev*archiv_sed[0], archiv_sed[2],xerr = archiv_sed[1],yerr = archiv_sed[3],color='r',ls ="None",fmt =".",label="Archived")

# 6 point SED of a flare (data based on fermi analysis)
#ax.errorbar(flare.Energy_Center*10**6, flare.Flux,yerr=flare.Flux_Error,xerr =[fl_err1*10**6,fl_err2*10**6]
#            ,color='b',fmt = ".",label="Flaring state")


# 6 point SED for quit state (data based on fermi analysis)
#ax.errorbar(quit.Energy_Center*10**6, quit.Flux,yerr=quit.Flux_Error,xerr =[qt_err1*10**6,qt_err2*10**6],color='g',fmt = ".",label="Quiescent state")


# Power law functions

#Electron energy distribution can be power law with exponential cut off or broken power law.

def PowerLawExpCutOff(alpha, gamma_cutOff, gamma,**kwargs):

    return N0 * (gamma ** (-alpha)) * np.exp(-(gamma / gamma_cutOff))

def BrokenPowerLaw(alpha_1,alpha_2,gamma_break,gamma,**kwargs):
    if gamma < gamma_break:
        return N0*gamma**(-alpha_1)
    else:
        return gamma_break**(alpha_2-alpha_1)*N0*gamma**(-alpha_2)

# Choose from two laws with booleans
def law_selection(cutOff_bool,broken_bool,alpha,alpha_1,alpha_2,gamma,gamma_cutOff,gamma_break):
    if (cutOff_bool ==1) and (broken_bool == 0):
        return PowerLawExpCutOff(alpha, gamma_cutOff, gamma)
    elif (cutOff_bool ==0) and (broken_bool == 1):
        return  BrokenPowerLaw(alpha_1,alpha_2,gamma_break,gamma)

# This function shows the energy which have many emitted photons(characteristic energy)
def char_energy(B, gamma):
    return (3 * e * h * B * (gamma ** 2)) / (2 * m * c)


# Bessel's function
def bessel_function(x):
    return (2.15 * (x) ** (1 / 3)) * (
                ((1 + 3.06 * (x)) ** (1 / 6)) * (1 + (0.884 * (x) ** (2 / 3)) + (0.471 * (x) ** (4 / 3)))) / (
                       1 + (1.64 * (x) ** (2 / 3)) + (0.974 * (x) ** (4 / 3))) * np.exp((-x))


# Quantity of emitted photons in per second with E energy(Differential spectrum dN/dE*dt)
def diff_spectrum(photon_eng, B, gamma):
    return ((np.sqrt(3) * (e ** 3) * B) / (2 * np.pi * restenergy * h * photon_eng)) * bessel_function(
       photon_eng / char_energy(B, gamma))


# spectrum for a group of electrons
def emission_spectrum(B, alpha,alpha_1,alpha_2,photon_eng, gamma_cutOff,gamma_break, cutOff_bool,broken_bool):
    return \
    quad(lambda j: (10 ** j) * (np.log(10)) * diff_spectrum(photon_eng, B, 10 ** j) * law_selection(cutOff_bool,broken_bool,alpha,alpha_1,alpha_2,10 ** j,gamma_cutOff,gamma_break),
         np.log10(gamma_min), np.log10(gamma_max))[0]



# Bolometric Synchrotron luminosity wich we get multiplying emission spectrum by square of photon energy(erg s^-1)
def luminosity(B, alpha,alpha_1,alpha_2,photon_eng, gamma_cutOff,gamma_break, cutOff_bool,broken_bool):
    return photon_eng**2*emission_spectrum(B, alpha,alpha_1,alpha_2,photon_eng, gamma_cutOff,gamma_break, cutOff_bool,broken_bool)



# Energy on a per square in a per second will be luminosity divided to surface F = L/4*pi*R^2.
# Considering Lorentz transformation formulas F should be multiplied by doppler factor ^4. F = delta^4*F'
# As an energy unit are being used erg(for converting use ev to erg ratio).
def flux_our_system(B, alpha,alpha_1,alpha_2,photon_eng, gamma_cutOff,gamma_break, cutOff_bool,broken_bool):
    return (doppler_factor**2)*1/(evtoerg*distance_surf)*luminosity(B, alpha,alpha_1,alpha_2,photon_eng*(1+red_shift)/doppler_factor, gamma_cutOff,gamma_break, cutOff_bool,broken_bool)





def synchrotron_plotter(B, alpha, alpha_1, alpha_2, gamma_cutOff, gamma_break, cutOff_bool, broken_bool):

    #
    energy_axis = np.logspace(-9, 4, num=100)
    # for  each point of energy_axis should be calculated synch_flux(as a Y axis)
    synchrotron_flux = np.array(
        [flux_our_system(B, alpha, alpha_1, alpha_2, i, gamma_cutOff, gamma_break, cutOff_bool, broken_bool) for i in
         energy_axis])

    ax.plot(energy_axis, synchrotron_flux, color="gray")
    ax.plot()

synchrotron_plotter(0.02, 2.697877, None, None,(187.518*10**9)/(restenergy), None, 1, 0)
ax.set_title("SED of PKS 1430-178", fontsize=14, fontweight='bold')
ax.set_xlabel('E [eV]')
ax.set_ylabel(r"$\nu F(\nu)\/ (erg\/cm^{-2} s^{-1})$")
ax.set_ylim(10**-18,2*10**-14)
ax.set_xlim(10**-6,10**13)
ax.set_xscale('log')
ax.set_yscale('log')
#ax.legend()
plt.savefig('SED.pdf')
#plt.show()

