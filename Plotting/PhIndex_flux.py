import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

param_table = pd.read_csv("param_table.csv")
fermi_data = pd.DataFrame(param_table)


fig,ax = plt.subplots(figsize =(16,9))
for i in range(35):
    if fermi_data.Type[i] == "fsrq":
        ax.scatter(fermi_data.flux[i],fermi_data.PL_index[i],c="r",marker='^',label = "FSRQ")
    elif  fermi_data.Type[i] == "bcu":
        ax.scatter(fermi_data.flux[i], fermi_data.PL_index[i],c= "g",marker='o',label = "BCU")
    elif fermi_data.Type[i] == "bll":
        ax.scatter(fermi_data.flux[i], fermi_data.PL_index[i],c= "b",marker="*",label = "BLL")

ax.set_xlabel(r"$Flux\/[ph\ \/cm^{-2} s^{-1}] $")
ax.set_ylabel("Photon index")


for i in range(35):
    ax.annotate(fermi_data.Assoc_name[i], xy=(fermi_data.flux[i], fermi_data.PL_index[i]),  xycoords='data',
       xytext=(0.01,6), textcoords='offset points')
plt.savefig("Flux-Index.pdf")
plt.show()



