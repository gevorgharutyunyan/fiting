import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

h = 6.5821 * 10 ** (-16)  # Plank's constant with line in Ev
hev = 4.135667662 * 10 ** -15  # Plank's constant  in Ev

archiv_data = pd.read_csv("PKS1430178.txt",header=None,delim_whitespace=True,error_bad_lines=False,warn_bad_lines=False)
archiv_sed = pd.DataFrame(archiv_data)

flare_data = pd.read_csv("flare_pks_1430_178_sed_data.csv",delim_whitespace=True,
                                                                             error_bad_lines=False,warn_bad_lines=False)
flare = pd.DataFrame(flare_data)
flare=flare[flare.Flux_Error !=0]

quit_data = pd.read_csv("pks_1430_178_sed_data.csv",delim_whitespace=True,
                                                                             error_bad_lines=False,warn_bad_lines=False)
quit = pd.DataFrame(quit_data)
quit=quit[quit.Flux_Error !=0]

fl_err1 = flare.Energy_Center-flare.Energy_Min
fl_err2 = flare.Energy_Max-flare.Energy_Center

qt_err1 = quit.Energy_Center-quit.Energy_Min
qt_err2 = quit.Energy_Max-quit.Energy_Center



fig, ax = plt.subplots(figsize=(14, 5))

# Archived data from asi.it
ax.errorbar(hev*archiv_sed[0], archiv_sed[2],xerr = archiv_sed[1],yerr = archiv_sed[3],color='r',ls ="None",fmt =".",label="Archived")

# 6 point SED of a flare (data based on fermi analysis)
ax.errorbar(flare.Energy_Center*10**6, flare.Flux,yerr=flare.Flux_Error,xerr =[fl_err1*10**6,fl_err2*10**6]
            ,color='b',fmt = ".",label="Flaring state")


# 6 point SED for quit state (data based on fermi analysis)
ax.errorbar(quit.Energy_Center*10**6, quit.Flux,yerr=quit.Flux_Error,xerr =[qt_err1*10**6,qt_err2*10**6],color='g',fmt = ".",label="Quiescent stat")


ax.set_title("SED of PKS 1430-178", fontsize=14, fontweight='bold')
ax.set_xlabel('E [eV]')
ax.set_ylabel(r"$\nu F(\nu)\/ (erg\/cm^{-2} s^{-1})$")
ax.set_ylim(10**-14,2*10**-10)
ax.set_xlim(10**-6,10**13)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.show()

