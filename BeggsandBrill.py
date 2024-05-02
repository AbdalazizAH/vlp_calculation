import math
import Fluid


def Flow_regime(Nfr, laml, L1, L2, L3, L4):
    """Function to Determine the Flow Regime by the Method of Beggs and Brill"""

    if ((laml < 0.01) and (Nfr < L1)) or ((laml >= 0.01) and (Nfr < L2)):
        flow_regime = 1

    if (laml >= 0.01) and (L2 < Nfr) and (Nfr <= L3):
        flow_regime = 2

    if (((0.01 <= laml) and (laml < 0.4)) and ((L3 < Nfr) and (Nfr < L1))) or (
        (laml >= 0.4) and (L3 < Nfr) and (Nfr <= L4)
    ):
        flow_regime = 3

    if ((laml < 0.4) and (Nfr >= L1)) or ((laml >= 0.4) and (Nfr > L4)):
        flow_regime = 4

    return flow_regime


def Liq_holdup(Nfr, Nvl, laml, angle, regime):
    """Function to Calculate Liquid Holdup for the Segregated, Intermittent and Distributed Regimes
    by the Method of Beggs and Brill"""
    if regime == 1:
        a_val = 0.98
        b_val = 0.4846
        c_val = 0.0868
        if angle >= 0:
            d_val = 0.011
            e_val = -3.768
            f_val = 3.539
            g_val = -1.614
        else:
            d_val = 4.7
            e_val = -0.3692
            f_val = 0.1244
            g_val = -0.5056

    if regime == 3:
        a_val = 0.845
        b_val = 0.5351
        c_val = 0.0173
        if angle >= 0:
            d_val = 2.96
            e_val = 0.305
            f_val = -0.4473
            g_val = 0.0978
        else:
            d_val = 4.7
            e_val = -0.3692
            f_val = 0.1244
            g_val = -0.5056

    if regime == 4:
        a_val = 1.065
        b_val = 0.5824
        c_val = 0.0609
        if angle >= 0:
            d_val = 1
            e_val = 0
            f_val = 0
            g_val = 0
        else:
            d_val = 4.7
            e_val = -0.3692
            f_val = 0.1244
            g_val = -0.5056

    corr = (1 - laml) * math.log(d_val * laml**e_val * Nvl**f_val * Nfr**g_val)
    if corr < 0:
        corr = 0

    psi = 1 + corr * (math.sin(1.8 * angle) - (math.sin(1.8 * angle)) ** 3 / 3)
    ylo = a_val * laml**b_val / Nfr**c_val
    if ylo < laml:
        ylo = laml

    yl_val = ylo * psi
    return yl_val


def Fric(Nre, eps):
    """Calculate Fanning Friction Factor using the Chen Equation"""
    try:
        math.log
        Temp = -4 * math.log10(
            (eps / 3.7065)
            - (5.0452 / Nre)
            * math.log10((eps**1.1098 / 2.8257) + (7.149 / Nre) ** 0.8981)
        )
    except Exception as inst:
        print(type(inst))
        print(inst.args)
        print(inst)

    return (1 / Temp) ** 2


def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """Calculate the Flowing Pressure Gradient by the Method of Beggs and Brill"""
    pi = math.pi
    Psep = 114.7
    Tsep = 50
    angle_rad = angle * pi / 180

    Z = Fluid.zfact((T + 460) / Fluid.Tc(gas_grav), P / Fluid.Pc(gas_grav))
    Wor = wtr_rate / oil_rate
    TDS = Fluid.salinity(wtr_grav)
    Pb = Fluid.Pbub(T, Tsep, Psep, gas_grav, oil_grav, Gor)
    Rso = Fluid.sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav)
    Rsw = Fluid.sol_gwr(P, T, TDS)
    Bo = Fluid.oil_fvf(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)
    Bw = Fluid.wtr_fvf(P, T, TDS)
    Bg = Fluid.gas_fvf(P, T, gas_grav)
    muo = Fluid.oil_visc(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)
    muw = Fluid.wtr_visc(P, T, TDS)
    mug = Fluid.gvisc(P, T + 460, Z, gas_grav)
    rhoo = Fluid.oil_dens(T, P, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)
    rhow = 62.368 * wtr_grav / Bw
    rhog = 2.699 * gas_grav * P / (T + 460) / Z
    sigo = Fluid.oil_tens(P, T, oil_grav)
    sigw = Fluid.wtr_tens(P, T)

    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)
    mul = (Bw * Wor * rhow) / (Bw * Wor * rhow + Bo * rhoo) * muw + (Bo * rhoo) / (
        Bw * Wor * rhow + Bo * rhoo
    ) * muo
    sigl = (Bw * Wor * rhow) / (Bw * Wor * rhow + Bo * rhoo) * sigw + (Bo * rhoo) / (
        Bw * Wor * rhow + Bo * rhoo
    ) * sigo

    qo = Bo * oil_rate / 15387
    qw = Bw * Wor * oil_rate / 15387
    ql = qo + qw
    qg = Bg * (Gor - Rso - Rsw * Wor) * oil_rate / 86400

    Axs = pi / 4 * (d / 12) ** 2
    usl = ql / Axs
    usg = qg / Axs
    um = usl + usg

    Nfr = um**2 / (d / 12) / 32.174
    Nvl = 1.938 * usl * (rhol / sigl) ** 0.25
    laml = usl / um
    lamg = 1 - laml
    L1 = 316 * laml**0.302
    L2 = 0.0009252 * laml**-2.4684
    L3 = 0.1 * laml**-1.4516
    L4 = 0.5 * laml**-6.738

    regime = Flow_regime(Nfr, laml, L1, L2, L3, L4)

    if regime == 2:
        a_val = (L3 - Nfr) / (L3 - L2)
        yl_seg = Liq_holdup(Nfr, Nvl, laml, angle_rad, 1)
        yl_int = Liq_holdup(Nfr, Nvl, laml, angle_rad, 3)
        yl_val = a_val * yl_seg + (1 - a_val) * yl_int
    else:
        yl_val = Liq_holdup(Nfr, Nvl, laml, angle_rad, regime)

    yg = 1 - yl_val

    rhom = rhol * laml + rhog * lamg
    mum = mul * laml + mug * lamg
    rhobar = rhol * yl_val + rhog * yg

    Nre = 1488 * rhom * um * (d / 12) / mum
    fn = Fric(Nre, 0.0006)
    x = laml / yl_val**2
    if 1 < x < 1.2:
        s = math.log(2.2 * x - 1.2)
    else:
        s = math.log(x) / (
            -0.0523
            + 3.182 * math.log(x)
            - 0.8725 * (math.log(x)) ** 2
            + 0.01853 * (math.log(x)) ** 4
        )

    ftp = fn * math.exp(s)

    Pgrad_pe = rhobar * math.sin(angle_rad) / 144
    Pgrad_f = 2 * ftp * rhom * um**2 / 32.17 / (d / 12) / 144
    Ek = um * usg * rhobar / 32.17 / P / 144

    return (Pgrad_pe + Pgrad_f) / (1 - Ek)
