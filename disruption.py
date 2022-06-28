import astropy.constants as cc
import numpy as np

msun_si = cc.M_sun.si.value
g_si = cc.G.si.value
c_si = cc.c.si.value
mass_conversion_constant = 0.075

def calculate_r_isco(dimensionless_spin,mass,type_):
    """
    Calculate BH isco radius

    :param dimensionless_spin: dimensionless spin magnitude
    :param mass: mass of the black hole in solar masses,
    :param type_: co_rotating or counter_rotating orbit
    :return: radius in km
    """
    a_tilde = dimensionless_spin * mass
    Z_1 = 1 + (1-a_tilde**2/mass**2)**(1/3) * ((1+a_tilde/mass)**(1/3) + (1-a_tilde/mass)**(1/3))
    Z_2 = (3*a_tilde**2/mass**2 + Z_1**2)**(1/2)

    if type_=='co_rotating':
        r_isco = mass * (3 + Z_2 - ((3-Z_1)*(3+Z_1+2*Z_2))**(1/2))
    elif type_=='counter_rotating':
        r_isco = mass * (3 + Z_2 + ((3-Z_1)*(3+Z_1+2*Z_2))**(1/2))
    r_isco_si = r_isco * msun_si * g_si / c_si ** 2
    return r_isco_si / 1000

def calculate_r_disruption(mbh, mns, rns):
    """
    Calculate disruption radius

    :param mbh:  mass of the black hole in solar masses,
    :param mns:  mass of the NS in solar masses,
    :param rns:  radius of NS in KM,
    :return: Disruption radius in km
    """
    return rns*(mbh/mns)**(1/3)

def nsbh_disruption(mbh, mns, rns, dimensionless_spin, type):
    """
    Calculate whether there is disruption or not

    :param mbh:  mass of the black hole in solar masses,
    :param mns:  mass of the NS in solar masses,
    :param rns:  radius of NS in KM,
    :param dimensionless_spin: dimensionless spin magnitude of BH
    :param type_: co_rotating or counter_rotating orbit
    :return: Boolean for disruption or not
    """
    rdis = calculate_r_disruption(mbh=mbh, mns=mns, rns=rns)
    risco = calculate_r_isco(dimensionless_spin=dimensionless_spin, mass=mbh, type_=type)
    return rdis >= risco

def rest_mass_from_gravitational_mass(gravitational_mass):
    """
    Convert from gravitational_mass to rest mass

    :param gravitational_mass: gravitational mass
    :return: rest mass
    """
    rest_mass = gravitational_mass + mass_conversion_constant * gravitational_mass ** 2
    return rest_mass


def gravitational_mass_from_rest_mass(rest_mass):
    """
    Convert from rest mass to gravitational mass

    :param rest_mass: rest/baryonic mass
    :return: gravitational mass
    """
    gravitational_mass = (-1. + np.sqrt(1 + 4. * mass_conversion_constant * rest_mass)) / (
                2. * mass_conversion_constant)
    return gravitational_mass


def total_gravitational_mass(gravitational_mass_1, gravitational_mass_2, ejecta_mass=0.05):
    """
    Calculate remnant mass of BNS merger product

    :param gravitational_mass_1: mass of primary ns in solar masses
    :param gravitational_mass_2: mass of secondary ns in solar masses
    :param ejecta_mass: mass ejected from system in solar masses
    :return: remnant mass
    """
    rest_mass_1 = rest_mass_from_gravitational_mass(gravitational_mass_1)
    rest_mass_2 = rest_mass_from_gravitational_mass(gravitational_mass_2)
    total_rest_mass = rest_mass_1 + rest_mass_2 - ejecta_mass
    total_gravitational_mass = gravitational_mass_from_rest_mass(total_rest_mass)
    return total_gravitational_mass

def bns_disruption(m1, m2, mtov, ejecta_mass=0.05):
    """

    :param m1: mass of primary ns in solar masses
    :param m2: mass of secondary ns in solar masses
    :param mtov: maximum allowed non-rotating mass, in solar masses
    :param ejecta_mass: ejecta in solar masses
    :return: boolean for whether BNS can launch jet.
    """
    remnant_mass = total_gravitational_mass(m1, m2, ejecta_mass=ejecta_mass)
    return remnant_mass >= 1.2 * mtov
