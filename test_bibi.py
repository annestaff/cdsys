from math import exp

import matplotlib.pyplot as plt
import numpy as np
from ddeint import ddeint

V0 = 1e4
T0 = 7e7


def model(Y, t, pv=210, delta_v=5, kl=0.4, beta=5e-7, gt=0.8, beta_prim=3e-8, rho=2.6, fi=0.33, delta_i=2,
          kn=2.5, ke=5e-5, pf=1e-5, delta_f=2, beta_cn=1, hc=1e4, tau_c=6, pc=1.2, delta_e=0.57, beta_bn=0.03, hb=1e4,
          tau_b=4, delta_p=0.5, pb=0.52, ps=12, delta_s=2, pl=4, delta_l=0.015, ks=0.8):
    viral_load, target_cells, infected_cells, ifn, virus_resistance, naive_CD8, effector, naive_B, plasma_cells, short_antibodies, long_antibodies = Y(t)

    viral_load_tau_c, naive_CD8_tau_c = Y(t - tau_c)[0], Y(t - tau_c)[5]
    viral_load_tau_b, naive_B_tau_b = Y(t - tau_b)[0], Y(t - tau_c)[7]

    VIRAL_LOAD = pv * infected_cells - delta_v * viral_load - ks * viral_load * short_antibodies - kl * viral_load * long_antibodies - beta * viral_load * target_cells

    TARGET_CELLS = gt * (target_cells + virus_resistance) * (1 - (target_cells + virus_resistance + infected_cells) / T0) - beta_prim * viral_load * target_cells + rho * virus_resistance - fi * ifn * target_cells

    INFECTED_CELLS = beta_prim * viral_load * target_cells - delta_i * infected_cells - kn * infected_cells * ifn - ke * infected_cells * effector

    IFN = pf * infected_cells - delta_f * ifn

    VIRUS_RESISTANCE = fi * ifn * target_cells - rho * virus_resistance

    NAIVE_CD8 = - beta_cn * (viral_load / (viral_load + hc)) * naive_CD8

    if t < tau_c:
        EFFECTOR = 0
    else:
        EFFECTOR = beta_cn * (viral_load_tau_c / (viral_load_tau_c + hc)) * naive_CD8_tau_c * exp(pc * tau_c) - delta_e * effector

    NAIVE_B = - beta_bn * (viral_load / (viral_load + hb)) * naive_B

    if t < tau_b:
        PLASMA_CELLS = 0
    else:
        PLASMA_CELLS = beta_bn * (viral_load_tau_b / (viral_load_tau_b + hb)) * naive_B_tau_b * exp(pb * tau_b) - delta_p * plasma_cells

    SHORT_ANTIBODIES = ps * plasma_cells - delta_s * short_antibodies

    LONG_ANTIBODIES = pl * plasma_cells - delta_l * long_antibodies

    return np.array(
        [VIRAL_LOAD, TARGET_CELLS, INFECTED_CELLS, IFN, VIRUS_RESISTANCE, NAIVE_CD8, EFFECTOR, NAIVE_B, PLASMA_CELLS,
         SHORT_ANTIBODIES, LONG_ANTIBODIES])


legend = ['V', 'T', 'I', 'F', 'R', 'Cn', 'E', 'Bn', 'P', 'As', 'Al']


def run_model(y0=np.array([V0, T0, 0, 0, 0, 100, 0, 100, 0, 0, 0]), t0=0, interval=12, graph_step=10000):
    tt = np.linspace(t0, interval, graph_step)

    y = ddeint(model, lambda t: y0, tt)
    return {"y0": y0, "tt": tt, "y": y}


def plot_model(y, tt=np.linspace(0, 20, 10000)):
    plt.close('all')
    plt.figure(1)
    plt.clf()
    plt.subplot(531)
    plt.plot(tt, y[:, 0], label=legend[0])
    plt.legend()
    plt.subplot(532)
    plt.plot(tt, y[:, 1], label=legend[1])
    plt.legend()
    plt.subplot(533)
    plt.plot(tt, y[:, 2], label=legend[2])
    plt.legend()
    plt.subplot(534)
    plt.plot(tt, y[:, 3], label=legend[3])
    plt.legend()
    plt.subplot(535)
    plt.plot(tt, y[:, 4], label=legend[4])
    plt.legend()
    plt.subplot(536)
    plt.plot(tt, y[:, 5], label=legend[5])
    plt.legend()
    plt.subplot(537)
    plt.plot(tt, y[:, 6], label=legend[6])
    plt.legend()
    plt.subplot(538)
    plt.plot(tt, y[:, 7], label=legend[7])
    plt.legend()
    plt.subplot(539)
    plt.plot(tt, y[:, 8], label=legend[8])
    plt.legend()
    plt.subplot(5310)
    plt.plot(tt, y[:, 9], label=legend[9])
    plt.legend()
    plt.subplot(5311)
    plt.plot(tt, y[:, 10], label=legend[10])
    plt.legend()
    plt.show()


def plot_model_2(y, y0=np.array([V0, T0, 0, 0, 0, 100, 0, 100, 0, 0, 0]), tt=np.linspace(0, 12, 10000)):
    for u in range(len(y0)):
        plt.plot(tt, y[:, u], label=legend[u])
        plt.legend()
        plt.show()


my_model = run_model()
plot_model_2(my_model["y"])
print("Done")
