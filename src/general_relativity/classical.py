import math

G = 6.67430e-11
c = 299792458.0


def deflection_angle_point_mass(M, b):
    """GR light deflection (radians) for a point mass in weak field: alpha = 4GM/(c^2 b)."""
    return 4 * G * M / (c * c * b)


def shapiro_delay_point_mass(M, r1, r2, b):
    """
    Shapiro time delay (seconds) for signal passing with impact parameter b near mass M.
    Uses standard logarithmic formula: dt = (2GM/c^3) ln[(r1 + r2 + D)/(r1 + r2 - D)] with D≈sqrt(r1^2 + r2^2 - 2 r1 r2 cos(theta)).
    For simplicity assume small-angle path with closest approach b and large r1,r2: dt ≈ (2GM/c^3) ln(4 r1 r2 / b^2).
    """
    return (2 * G * M / (c ** 3)) * math.log(4 * r1 * r2 / (b * b))


def gravitational_redshift_potential(delta_phi):
    """Gravitational redshift z ≈ delta_phi/c^2 for small potentials (frequency shift)."""
    return delta_phi / (c * c)


def perihelion_precession(a, e, M):
    """Per orbit perihelion advance (radians): delta_omega = 6*pi*GM/(c^2 a (1-e^2))."""
    return 6 * math.pi * G * M / (c * c * a * (1 - e * e))


def lense_thirring_precession(J, r):
    """Frame-drag angular velocity omega = 2 G J / (c^2 r^3)."""
    return 2 * G * J / (c * c * r ** 3)


def newtonian_acceleration(M, r_vec):
    """Newtonian acceleration vector: a = - G M r / |r|^3."""
    rx, ry, rz = r_vec
    r3 = (rx * rx + ry * ry + rz * rz) ** 1.5
    if r3 == 0:
        return (0.0, 0.0, 0.0)
    scale = -G * M / r3
    return (scale * rx, scale * ry, scale * rz)


def einstein_radius_point_mass(M, D_l, D_s, D_ls):
    """Einstein angle (radians): theta_E = sqrt(4GM D_ls / (c^2 D_l D_s))."""
    return math.sqrt(4 * G * M * D_ls / (c * c * D_l * D_s))


def gw_phase_speed():
    """GW phase/group speed in GR vacuum: c."""
    return c


def binary_quadrupole_power(M1, M2, a):
    """Quadrupole power radiated by a circular binary: P = (32/5) G^4 (M1 M2)^2 (M1+M2) / (c^5 a^5)."""
    return (32.0 / 5.0) * (G ** 4) * ((M1 * M2) ** 2) * (M1 + M2) / (c ** 5 * a ** 5)


def circular_orbit_energy(M, a):
    """Specific orbital energy (per unit mass) for circular orbit in Newtonian limit: E = - GM/(2a)."""
    return -G * M / (2 * a)


def binary_period_derivative(M1: float, M2: float, P_b: float, e: float) -> float:
    """Peters (1964) orbital period derivative for an eccentric binary (GR).

    dP_b/dt = -(192π/5) (G/c³)^(5/3) (2π/P_b)^(5/3) M1 M2 / (M1+M2)^(1/3) f(e)

    where f(e) = (1 + 73e²/24 + 37e⁴/96) / (1-e²)^(7/2).
    """
    import math
    e2 = e * e
    f_e = (1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2) / (1.0 - e2) ** 3.5
    M_tot = M1 + M2
    return (
        -(192.0 * math.pi / 5.0)
        * (G / c ** 3) ** (5.0 / 3.0)
        * (2.0 * math.pi / P_b) ** (5.0 / 3.0)
        * M1 * M2 / M_tot ** (1.0 / 3.0)
        * f_e
    )


def gr_surface_redshift(M_star: float, R_star: float) -> float:
    """GR (Schwarzschild) surface gravitational redshift.

    z_GR = 1 / sqrt(1 − 2GM/c²R) − 1

    Parameters
    ----------
    M_star : float   Total stellar mass [kg].
    R_star : float   Stellar radius [m].

    Returns
    -------
    float   Dimensionless surface redshift.
    """
    compactness = 2.0 * G * M_star / (c * c * R_star)
    return 1.0 / (1.0 - compactness) ** 0.5 - 1.0


# Leaver (1985) l=2, n=0 Schwarzschild QNM coefficients (dimensionless, × GM/c³)
_QNM_OMEGA_R = 0.37367   # Re(ω) × GM/c³  — oscillation
_QNM_OMEGA_I = 0.08896   # Im(ω) × GM/c³  — inverse decay (half-rate)


def gr_qnm_frequency(M_bh: float) -> float:
    """GR Schwarzschild QNM ring frequency for the fundamental (l=2, n=0) mode.

    f_QNM = ω_R / (2π)   where  ω_R = 0.37367 × c³/(GM)

    From Leaver (1985) exact continued-fraction solutions; tabulated in
    Berti, Cardoso & Starinets (2009), Table I.

    Parameters
    ----------
    M_bh : float   Black-hole mass [kg].

    Returns
    -------
    float   Ring-down oscillation frequency [Hz].
            Scales as 12.1 kHz × (M_sun / M_bh).
    """
    import math
    return _QNM_OMEGA_R * c ** 3 / (2.0 * math.pi * G * M_bh)


def gr_qnm_damping_time(M_bh: float) -> float:
    """GR Schwarzschild QNM e-folding damping time for the fundamental (l=2, n=0) mode.

    τ_QNM = 1 / (2 ω_I)   where  ω_I = 0.08896 × c³/(GM)

    Parameters
    ----------
    M_bh : float   Black-hole mass [kg].

    Returns
    -------
    float   Damping e-folding time [s].  Scales as 27.8 μs × (M_bh / M_sun).
    """
    return G * M_bh / (2.0 * _QNM_OMEGA_I * c ** 3)


def gr_photon_sphere_radius(M_bh: float) -> float:
    """GR Schwarzschild photon-sphere radius: r_ps = 3 GM/c².

    Parameters
    ----------
    M_bh : float   Black-hole mass [kg].

    Returns
    -------
    float   Photon-sphere radius [m].
    """
    return 3.0 * G * M_bh / (c * c)
