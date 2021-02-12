import numpy as np
from create_rayleigh_curves import RayleighTropen

def aeq(T):
    return 1 / np.exp( 24844 / T**2 - 76.248 / T + 0.05261 )

def get_dDdeq(T, dDv, dDr):
    """
    INPUT:
    T   : temperature (in degree Celsius)
    dDv : dD of vapor
    dDr : dD of rain

    RETURN:
    dDdeq : disequilibrium between vapor and rain
    """
    
    dDeq  = aeq(T + 273.15) * (dDr + 1000) - 1000
    dDdeq = dDv - dDeq

    return dDdeq


def calc_rayleigh(T0, rh0, dD0, f = 1, T_g = -30):
    """
    INPUT:
    T0  : starting temperature, given degree Celsius
    rh0 : starting relative humidity, given in percent
    dD0 : starting dD, given in permil
    f   : factor for adjusting alpha
    T_g : boundary between liquid and frozen condensation

    RETURN:
    h2o_rl : H2O along Rayleigh process
    dD_rl  : dD along Rayleigh process
    d18O   : d18O along Rayleigh process
    t_rl   : temperature along Rayleigh process
    p_rl   : pressure along Rayleigh process
    """
    return RayleighTropen(T0, rh0, dD0, f, T_g)

def calc_mixing(p0, p1, l_log = True):
    """
    INPUT:
    p0 : dry mixing member
         list: [H2O, dD]
    p1 : moist mixing member
         list: [H2O, dD]
    
    RETURN:
    q  : H2O along mixing curve
    dD : dD along mixing curve
    """
    q0, dD0 = p0
    q1, dD1 = p1

    if l_log:
        q = np.exp(np.linspace(np.log(q0), np.log(q1)))
    else:
        q = np.linspace(q0, q1)
    dDf = (dD1 * q1 - dD0 * q0) / (q1 - q0)
    dD = q0 * (dD0 - dDf) * (1/q) + dDf

    return q, dD


def plot_rayleigh(ax, rl_list = None, lwf = 1, t = None, rh = None, dD0 = None, c = 'dimgray'):
    """
    INPUT:
    ax: plot axis where to plot the Rayleigh curve
    rl_list = [ rl_line_1, rl_line_2, ... ] with rl_line_i = [T0 in Celsius, RH0, dD0]

    RETURN:
    cs_rl : line plot instance
    """
    if rl_list is None:
        if t is None and rh is None and dD0 is None:
            rl_list = [ [30, 90, -80] ]
        else:
            rl_list = [ [t, rh, dD0] ]
    
    for rl in rl_list:
        h2o_rl, dD_rl, d18_rl, t_rl, p_rl = calc_rayleigh(rl[0], rl[1], rl[2], T_g = -10)
        cs_rl, = ax.plot(h2o_rl, dD_rl, c = c, ls = (0,(1,1)), lw = 4.5*lwf,
                         #ls = ':', lw = 5,
                                             label = 'Rayleigh')
    return cs_rl

def plot_superrayleigh(ax, srl_list = None, lwf = 1, c = 'gray', lxl = None):
    """
    srl_list = [ srl_line_1, srl_line_2, ... ] with srl_line_i = [T0 in Celsius, RH0, dD0, alpha_fac]
    """
    if srl_list is None:
        srl_list = [ [13, 80,  -175, 1.5],
             [22, 100, -110, 1.3] ]

    for sx, srl in enumerate(srl_list):
        if lxl == None:
            if sx == 0:
                lx = 100
            else:
                lx = 130
        else:
            lx = lxl[sx]

        h2o_sr, dD_sr, _, _, _ = calc_rayleigh(srl[0], srl[1], srl[2], f = srl[3], T_g = -10)
        cs_srl, = ax.plot(h2o_sr[:lx], dD_sr[:lx], c = c, ls = (0, (3,1,1,1,1,1)), lw = 4.5*lwf, 
                          #ls = '--', lw = 5, 
                         label = 'Super-Rayleigh')
    return cs_srl

def plot_mixing(ax, mx_list = None, lwf = 1, c = 'gray'):
    """
    mx_list = [ mx_line_1, mx_line_2, ... ] with mx_line_i = [ [q1,dD1], [q2,dD2] ]
    """
    if mx_list is None:
        mx_list = [
             # [[6e2, -700], [3e4, -120]],
              [[5e1, -700], [1.53e4, -150]],
             # [[5e1, -700], [1.53e4, -120]],
              [[1e3, -450], [2.2e4, -120]],
            #  [[1e4, -250], [2.6e4, -120]],
              ]

    for mx in mx_list:
        h2o_mx, dD_mx = calc_mixing(p0 = mx[0], p1 = mx[1])
        cs_mx, = ax.plot(h2o_mx, dD_mx, c = c, ls = (0,(3,1)), lw = 4*lwf,
                         #(0, (3,1,1,1,1,1)), lw = 5,
                                             label = 'mixing')
    return cs_mx
