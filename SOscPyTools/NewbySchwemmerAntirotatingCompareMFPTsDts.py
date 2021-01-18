import matplotlib

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})


from matplotlib import pyplot as plt
from matplotlib import cm

from SOscFilePy import *
from pathlib import Path
from SimConfigPy import *

if __name__ == "__main__":

    simConfig = SimConfig()
    simConfig.load_from_file(Path("../configs/config_cluster.json"))
    print("Input: " + str(simConfig.Paths["In"].resolve()))
    print("Output: " + str(simConfig.Paths["Out"].resolve()))
    IsoSet = IsochroneSet(simConfig.Paths["In"].resolve(), \
            simConfig.ModelName)
    mfptSet = MFPTSet(simConfig.Paths["Out"].resolve(), \
            simConfig.ModelName)

    latex_textwidth_in = 6.50127
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(w=latex_textwidth_in, h=0.5*latex_textwidth_in)

    (pset, rhos, phis) = IsoSet.__get_Isochrone__("Isochrone6")
    I = mfptSet.__get_Isochrone__("Isochrone6")
    ax1.plot(I.Rho_init, I.MFPT, '.', color="blue", \
            label=r"dt = $" + str(simConfig.Simulation["dt"]) + \
            "$, N = $" + str(simConfig.Simulation["Ensemble Size"]) + "$")
    ax1.errorbar(I.Rho_init, I.MFPT, yerr=np.sqrt(I.VarFPT),\
            ls=" ", mec="blue")

    
    del IsoSet
    del mfptSet
    
    simConfig.load_from_file(Path("../configs/config_cluster_high_noise.json"))
    print("Input: " + str(simConfig.Paths["In"].resolve()))
    print("Output: " + str(simConfig.Paths["Out"].resolve()))
    IsoSet = IsochroneSet(simConfig.Paths["In"].resolve(), \
            simConfig.ModelName)
    mfptSet = MFPTSet(simConfig.Paths["Out"].resolve(), \
            simConfig.ModelName)

    (pset, rhos, phis) = IsoSet.__get_Isochrone__("Isochrone8")
    I = mfptSet.__get_Isochrone__("Isochrone0")
    ax1.plot(I.Rho_init, I.MFPT, '.', color="orange", \
            label=r"dt = $" + str(simConfig.Simulation["dt"]) + \
            "$, N = $" + str(simConfig.Simulation["Ensemble Size"]) + "$")
    ax1.errorbar(I.Rho_init, I.MFPT, yerr=np.sqrt(I.VarFPT),\
            ls=" ", mec="orange")
    
    del IsoSet
    del mfptSet 

    simConfig = SimConfig()
    simConfig.load_from_file(Path("../configs/config_cluster_high_noise2.json"))
    print("Input: " + str(simConfig.Paths["In"].resolve()))
    print("Output: " + str(simConfig.Paths["Out"].resolve()))
    IsoSet = IsochroneSet(simConfig.Paths["In"].resolve(), \
            simConfig.ModelName)
    mfptSet = MFPTSet(simConfig.Paths["Out"].resolve(), \
            simConfig.ModelName)

    (pset, rhos, phis) = IsoSet.__get_Isochrone__("Isochrone8")
    I = mfptSet.__get_Isochrone__("Isochrone0")
    ax1.plot(I.Rho_init, I.MFPT, '.', color="red", \
            label=r"dt = $" + str(simConfig.Simulation["dt"]) + \
            "$, N = $" + str(simConfig.Simulation["Ensemble Size"]) + "$")
    ax1.errorbar(I.Rho_init, I.MFPT, yerr=np.sqrt(I.VarFPT),\
            ls=" ", mec="red")

    ax2.plot(rhos*np.cos(phis), rhos*np.sin(phis), color="black")


    phis = np.linspace(0.0, 2*np.pi, num=100)
    rm = np.amin(rhos)*np.ones(phis.shape)
    rp = np.amax(rhos)*np.ones(phis.shape)
    ax2.plot(rm*np.cos(phis), rm*np.sin(phis), ":k")
    ax2.plot(rp*np.cos(phis), rp*np.sin(phis), "--k")


    ax2.axis("equal")
    ax2.set_title("Phasenraum")
    ax1.legend()
    ax1.set_xlabel(r"$\rho$")
    ax1.set_ylabel(r"MFPT")
    fig.suptitle(r"$D = " + str(pset["D"][0]) + "$")

    del IsoSet
    del mfptSet
    del simConfig

    plt.savefig("NewbySchwemmerCompareMFPTsDts.pgf")
    plt.savefig("NewbySchwemmerCompareMFPTsDts.pdf")


