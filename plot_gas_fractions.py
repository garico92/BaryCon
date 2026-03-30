import numpy as np
import matplotlib.pyplot as plt


#cosmology assumed for plotting purposes 

hubble = 0.6766 # Planck18 hubble constant 
f_uni =  0.04897/(0.04897+0.26067) #Planck18 universal baryon fraction, with respect to the cold matter (no massive neutrinos included) 

## X-rays datasets ##
###############################################################################################################
data_dir = "data/"

#WARNING: Cosmology dependence of the measurements to be taken into account

# Popesso 2024 bestfit relations, eROSITA + CHEXMATE, https://arxiv.org/pdf/2411.16555
h0_popesso = 0.6774
M200_ih = np.logspace(12.5,15)
M200 = M200_ih / h0_popesso * (h0_popesso/hubble)
fgas_200 = (2.09) * 1e-6 * (M200)**(0.33) * h0_popesso**(3/2)/(hubble**(3./2))
fgas_200_plus = (2.09+0.14) * 1e-6 * (M200)**(0.33 + 0.02) * h0_popesso**(3/2)/(hubble**(3./2))
fgas_200_minus = (2.09-0.14) * 1e-6 * (M200)**(0.33 - 0.02) * h0_popesso**(3/2)/(hubble**(3./2))

M500_ih = np.logspace(12.5,15)
M500 = M500_ih/h0_popesso * (h0_popesso/hubble)
fgas_500 = (2.23) * 1e-7 * (M500)**(0.39) * h0_popesso**(3/2)/(hubble**(3./2))
fgas_500_plus = (2.23+0.18) * 1e-7 * (M500)**(0.39 + 0.02) * h0_popesso**(3/2)/(hubble**(3./2))
fgas_500_minus = (2.23-0.18) * 1e-7 * (M500)**(0.39 - 0.02) * h0_popesso**(3/2)/(hubble**(3./2))

##### Bulbul et al. 2024 
data = np.loadtxt(data_dir+"Bulbul_2024.dat",skiprows=2)
h0_bulbul = 0.6715 #Hubble factor assumed in the paper
m500_B = data[:,0]*1e13*h0_bulbul 
m500_Bm = data[:,1]*1e13*h0_bulbul 
m500_Bp = data[:,2]*1e13*h0_bulbul 
f_gas_B = data[:,3] * h0_bulbul**(3/2)/(hubble**(3/2))
f_gas_Bm = data[:,4] * h0_bulbul**(3/2)/(hubble**(3/2))
f_gas_Bp = data[:,5] * h0_bulbul**(3/2)/(hubble**(3/2))

mask = f_gas_B >= 0

m500_B = m500_B[mask]
m500_Bm = m500_Bm[mask]
m500_Bp = m500_Bp[mask]
f_gas_B = f_gas_B[mask]
f_gas_Bm = f_gas_Bm[mask]
f_gas_Bp = f_gas_Bp[mask]

###### Akino regression  (Akino et al. 2022 https://arxiv.org/pdf/2111.10080, Appendix E)

hubble_Akino = 0.7
def Akino_regression(m500, alpha, beta):
    m = alpha + beta*np.log(m500/1e14) 
    M = np.exp(m)*1e12
    return M

M500_Akino = np.logspace(13, 15, 20)

M_gas_A = Akino_regression(M500_Akino, 2.02, 1.28)
M_gas_m1 = Akino_regression(M500_Akino, 2.02-0.10, 1.28-0.09)
M_gas_p1 = Akino_regression(M500_Akino, 2.02+0.08, 1.28+0.12)
M_gas_m2 = Akino_regression(M500_Akino, 2.02-0.10, 1.28+0.09)
M_gas_p2 = Akino_regression(M500_Akino, 2.02+0.08, 1.28-0.12)

f_gas_Akino = M_gas_A/M500_Akino
f_gas_Akino_m1 = M_gas_m1/M500_Akino
f_gas_Akino_p1 = M_gas_p1/M500_Akino
f_gas_Akino_m2 = M_gas_m2/M500_Akino
f_gas_Akino_p2 = M_gas_p2/M500_Akino
f_gas_Akino_m = np.min([f_gas_Akino_m1,f_gas_Akino_m2],axis=0)
f_gas_Akino_p = np.max([f_gas_Akino_p1,f_gas_Akino_p2],axis=0)

M500_Akino_ih = M500_Akino * hubble_Akino

############# Grandis et al. 2024 (https://arxiv.org/pdf/2309.02920)
data = np.genfromtxt(data_dir+"Grandis_2024.dat", comments="#",dtype=None, encoding=None)
h_Grandis = 0.67
# Columns: survey, component, Mpiv, fraction, frac_err, slope, slope_err
survey = data['f0']
component = data['f1']
M500_Grandis= data['f2'] * h_Grandis
fractions = data['f3']
err_frac = data['f4']
# Gas (ICM)
gas_mask = component == "ICM"
M_500c_gGrandis = M500_Grandis[gas_mask]
f_gas_Grandis = fractions[gas_mask]
err_f_gas_Grandis = err_frac[gas_mask]

############# Mantz et al. 2022 (https://arxiv.org/pdf/2111.09343)
data = np.genfromtxt(data_dir + "Mantz_2022.dat", dtype=None,
    names=True, encoding="utf-8", skip_header=16)
h_Mantz = 0.7
Om_Mantz = 0.3
z         = data["z"]
Mgas      = data["Mgas"]     * 1e13 * h_Mantz  # [Msun/h]
Mgas_err  = data["Mgas_err"] * 1e13 * h_Mantz 
Mtot      = data["Mtot"]     * 1e14  * h_Mantz  # [Msun/h]
Mtot_err  = data["Mtot_err"] * 1e14  * h_Mantz 
M2500     = data["M2500"]    * 1e14  * h_Mantz  # [Msun/h]
M2500_err = data["M2500_err"]* 1e14 * h_Mantz 
f_gas     = data["fgas"]
f_gas_err = data["fgas_err"]

###############Chiu et al. 2018 (https://academic.oup.com/mnras/article/478/3/3072/4996803)
h0_chiu = 0.68
data = np.genfromtxt(data_dir+"Chiu_2018.dat", dtype=None, names=True,encoding='utf-8', skip_header=17)#
M500_Chiu  = data['M500'] * 1e14 * h0_chiu
Mgas_Chiu  = data['MICM'] * 1e13 * h0_chiu
err_M500_Chiu = data['M500_err'] * 1e14 * h0_chiu
err_Mgas_Chiu = data['MICM_err'] * 1e13 * h0_chiu
f_gas_Chiu = Mgas_Chiu/M500_Chiu
err_f_gas_Chiu = np.sqrt((err_Mgas_Chiu/Mgas_Chiu)**2+(err_M500_Chiu/M500_Chiu)**2)*f_gas_Chiu
################################################################
##### Gonzalez (https://arxiv.org/pdf/1309.3565.pdf)
data = np.loadtxt(data_dir+"Gonzalez_2013.dat",skiprows=1)
h0_gonzalez = 0.6715 #Hubble factor assumed in the paper
f_gas = data[:,1] * h0_gonzalez**(3/2)/(hubble**(3/2))
erf_gas = data[:,2] * h0_gonzalez**(3/2)/(hubble**(3/2))
m500 = data[:,7]*1e14*h0_gonzalez 
erm500 = data[:,8]*1e14*h0_gonzalez
r500 = data[:,9]*h0_gonzalez
err500 = data[:,10]*h0_gonzalez
####### Vikhlinin (https://iopscience.iop.org/article/10.1086/500288/pdf)#WARNING: arXiv version (https://arxiv.org/pdf/astro-ph/0507092.pdf) is different!
data = np.loadtxt(data_dir+"Vikhlinin_2006.dat",skiprows=2)
h0_vik = 0.72 #Hubble factor assumed in the paper
f_gas_V = data[:,14] * h0_vik**(3/2)/(hubble**(3/2))
erf_gas_V = data[:,15] * h0_vik**(3/2)/(hubble**(3/2))
m500_V = data[:,10]*1e14 * h0_vik 
erm500_V = data[:,11]*1e14 * h0_vik
r500_V = data[:,0] * h0_vik
err500_V = data[:,1] * h0_vik
####### Giodini 2009 (https://arxiv.org/pdf/0904.0448.pdf)
data = np.loadtxt(data_dir+"Giodini_2009.dat",skiprows=2)
h0_gio = 0.72 #Hubble factor assumed in the paper ‹WANING: CHECK THIS
f_star_G = data[:,1]
erf_star_G = data[:,2]
m500_G = data[:,0] * h0_gio
f_gas_G = data[:,3] * h0_gio**(3/2)/(hubble**(3/2))
erf_gas_G = data[:,4] * h0_gio**(3/2)/(hubble**(3/2))
f_stargas_G = data[:,5]
erf_stargas_G = data[:,6]
####### Arnaud 2017 (https://www.aanda.org/articles/aa/pdf/2007/42/aa8541-07.pdf)
data = np.loadtxt(data_dir+"Arnaud_2007.dat",skiprows=2)
h0_arnaud = 0.7
m500_A = data[:,4] *h0_arnaud*1e14
m500_Ap = data[:,5] *h0_arnaud*1e14
m500_Am = -data[:,6] *h0_arnaud*1e14
f_gas_A = data[:,13] * h0_arnaud**(3/2)/(hubble**(3/2))
f_gas_Ap = data[:,14] * h0_arnaud**(3/2)/(hubble**(3/2))
f_gas_Am = -data[:,15] * h0_arnaud**(3/2)/(hubble**(3/2))
##########################################
#Mulroy et al. 2019 (https://arxiv.org/pdf/1901.11276)
h0_M = 0.7
Om_M = 0.3
data = np.genfromtxt(data_dir+"Mulroy_2019.dat", dtype=None, names=True, encoding='utf-8', skip_header=26)#
M_M         = data["M_post"]     * 1e14  * h0_M  # [M_sun/h]
M_M_ep      = data["M_post_err_plus"]  * 1e14 * h0_M
M_M_em      = data["M_post_err_minus"]  * 1e14 * h0_M
Mgas_M      = data["Mgas"]     * 1e14 * h0_M  # [M_sun/h]
Mgas_ep_M   = data["Mgas_err_plus"]  * 1e14 * h0_M
Mgas_em_M   = data["Mgas_err_minus"]  * 1e14 * h0_M
f_gas_M     = Mgas_M / M_M
f_gas_M_ep  = f_gas_M * np.sqrt((Mgas_ep_M / Mgas_M)**2 + (M_M_em / M_M)**2)   # upper
f_gas_M_em  = f_gas_M * np.sqrt((Mgas_em_M / Mgas_M)**2 + (M_M_ep / M_M)**2)   # lower

####### Sun 2008 (https://arxiv.org/pdf/0805.2320.pdf)
data = np.loadtxt(data_dir+"Sun_2008.dat",skiprows=2)
h0_sun = 0.73
m500_S = data[:,1]*h0_sun*1e13
m500_Sp = data[:,2]*h0_sun*1e13
m500_Sm = -data[:,3]*h0_sun*1e13
f_gas_S = data[:,4]* h0_sun**(3/2)/(hubble**(3/2))
f_gas_Sp = data[:,5]* h0_sun**(3/2)/(hubble**(3/2))
f_gas_Sm = -data[:,6]* h0_sun**(3/2)/(hubble**(3/2))
##########################################
########Ettori et al. 2011  https://arxiv.org/pdf/1009.3266
h0_ettori = 0.7
data = np.loadtxt(data_dir+"Ettori_2011_masses_method2.dat",skiprows=2) #M_200c
m200_ettori = data[:,6] * h0_ettori*1e14 
m200_ettori_p = data[:,7]*  h0_ettori*1e14
m200_ettori_m = -data[:,8]* h0_ettori*1e14
c200_ettori = data[:,3]

data = np.loadtxt(data_dir+"Ettori_2011.dat",skiprows=2) #f_500, method 2 
f_gas_ettori = data[:,-2] * h0_ettori**(3/2)/(hubble**(3/2)) 
f_gas_err_ettori = data[:,-1]* h0_ettori**(3/2)/(hubble**(3/2))


########Zhang et al. 2011 (https://arxiv.org/pdf/1109.0390)
data = np.genfromtxt(data_dir+"Zhang_2011.dat", dtype=None, names=True, encoding='utf-8',skip_header=1)
h0_Z = 0.7 #Hubble factor assumed in the paper
Om_Z= 0.3 #Omega matter assumed in the paper

M500_Z = data['M500'] * h0_Z*1e14
M500_err_Z = data['M500_err'] * h0_Z*1e14
f_gas_Z = data['fgas'] * h0_Z**(3/2)/(hubble**(3/2))
f_gas_err_Z = data['fgas_err'] * h0_Z**(3/2)/(hubble**(3/2))

def virial_mass_converter(delta1, mass1, conc1, delta2, conc2):
    """
    functions to compute halo mass considering an overdensity delta1 
    starting from a different overdensity delta2 and assuming a NFW halo profile
    e.g. how to go from Delta_200 to Delta_500
    see e.g. Ettori et al 2011
    """
    mass2 = delta2*conc2**3/(delta1*conc1**3)*mass1
    return mass2


def virial_concentration_converter(delta1, conc1, delta2):
    """
    functions to compute halo concentration
    changing value of overdensity and assuming a NFW halo profile
    e.g. how to go from Delta_200 to Delta_500
    see e.g. Ettori et al 2011
    """
    from scipy import optimize

    def function(conc2):
        return (conc1/conc2)**3 * (np.log(1+conc2)-conc2/(1+conc2)) / (
            np.log(1+conc1)-conc1/(1+conc1)) - delta2/delta1
    conc2 = optimize.fsolve(function, conc1)
    return conc2

c500_ettori = virial_concentration_converter(200,c200_ettori,500)
m500_ettori   = virial_mass_converter(200,m200_ettori,c200_ettori,500,c500_ettori)
m500_ettori_p =  virial_mass_converter(200,m200_ettori+m200_ettori_p,c200_ettori,500,c500_ettori) - m500_ettori
m500_ettori_m = m500_ettori-virial_mass_converter(200,m200_ettori-m200_ettori_m,c200_ettori,500,c500_ettori)


########Wicker et al 2022,  https://arxiv.org/pdf/2204.12823.pdf
h0_wicker = 0.7
data = np.loadtxt(data_dir+"Wicker_2022.dat",skiprows=2)#
m500_wicker = data[:,4] * h0_wicker*1e14
m500_wicker_p = data[:,5]*  h0_wicker*1e14
m500_wicker_m = data[:,6]* h0_wicker*1e14

f_gas_wicker = data[:,7]* h0_wicker**(3/2)/(hubble**(3/2))
f_gas_wicker_p = data[:,8]* h0_wicker**(3/2)/(hubble**(3/2))
f_gas_wicker_m = data[:,9]* h0_wicker**(3/2)/(hubble**(3/2))
##################################################################

######################## plotting #############################


fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(10,8), constrained_layout=True)
ax.axhline(f_uni,  c='k', label=r'$\Omega_{\rm b}/\Omega_{\rm m}$ Planck18', alpha=0.7)
c_pop = 'orange'
size=18
pad = 8

ax.errorbar(m500_V,f_gas_V,yerr=erf_gas_V,xerr=erm500_V,marker="o",c="red",ls="",label="Vikhlinin et al. 2006",alpha=0.9)
ax.errorbar(m500_S,f_gas_S,yerr=(f_gas_Sm,f_gas_Sp),xerr=(m500_Sm,m500_Sp),marker="o",c="g",ls="",label="Sun et al. 2008",alpha=0.9)
ax.errorbar(m500_G,f_gas_G,yerr=erf_gas_G, marker="o",c="blue",ls="",label="Giodini et al. 2009",alpha=0.9)
ax.errorbar(m500_ettori,f_gas_ettori,yerr=f_gas_err_ettori,xerr=(m500_ettori_m,m500_ettori_p), marker="*",c="purple",ls="",label="Ettori et al. 2011",alpha=0.9)
ax.errorbar(M500_Z, f_gas_Z, xerr=M500_err_Z, yerr=f_gas_err_Z, marker="o",c="cyan",ls="",label="Zhang et al. 2011",alpha=0.9)
ax.errorbar(m500,f_gas,yerr=erf_gas,xerr=erm500,marker="o",c="indigo",ls="",label="Gonzalez et al. 2013",alpha=0.9)
ax.errorbar(m500_A,f_gas_A,yerr=(f_gas_Am,f_gas_Ap),xerr=(m500_Am,m500_Ap),marker="o",c="m",ls="",label="Arnaud et al. 2017",alpha=0.9)
ax.errorbar(M500_Chiu,f_gas_Chiu,yerr=err_f_gas_Chiu,xerr=err_M500_Chiu, marker="o",c="violet",ls="",label="Chiu et al. 2018",alpha=0.9)
ax.errorbar(M_M,f_gas_M,yerr=(f_gas_M_em,f_gas_M_ep),xerr=(M_M_em,M_M_ep), marker="o",c="maroon",ls="",label="Mulroy et al. 2019",alpha=0.9)
ax.errorbar(m500_wicker,f_gas_wicker,yerr=(f_gas_wicker_m,f_gas_wicker_p),xerr=(m500_wicker_m,m500_wicker_p),marker="o",c="olive",ls="",label="Wicker et al. 2022",alpha=0.9)

ax.plot(m500_B,f_gas_B,marker="o",c="orange",ls="",alpha=0.05)
ax.plot([],[], marker="o",c="orange",ls="",label="Bulbul et al. 2024",alpha=0.5)
ax.semilogx(M500_ih, fgas_500, color=c_pop)
ax.fill_between(M500_ih,fgas_500_minus,fgas_500_plus,alpha=0.2, color=c_pop, label="Popesso et al. 2024")
ax.fill_between(M500_Akino_ih, f_gas_Akino_m, f_gas_Akino_p, color='brown', alpha=0.2, label="Akino et al. 2022")
ax.plot(M500_Akino_ih, f_gas_Akino, color='brown', alpha=0.9)
ax.errorbar(M_500c_gGrandis,f_gas_Grandis,yerr=err_f_gas_Grandis,marker="o",c="navy",ls="",label="Grandis et al. 2024",alpha=0.9)


ax.set_xlim(10**(12.5),2e15)
ax.set_ylim(0.,0.25)
ax.set_ylabel(r"$M_{\rm gas}/M_{\rm 500}$",size=size, rotation=90,labelpad=pad)
ax.set_xlabel(r"$M_{\rm 500} \, [h^{-1}\rm M_{\odot}]$",size=size, rotation=0,labelpad=pad)

ax.legend( prop={"size":14})#, fancybox=False, framealpha=0.,loc=‘center left’)

ax.tick_params(axis='both', which='major', pad=pad, labelsize=size)

plt.show()

