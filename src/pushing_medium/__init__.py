"""Pushing-Medium unified field theory implementation."""

from .core import (
    # Constants
    G, c, mu_0, epsilon_0, e, m_e,
    
    # Index field functions
    index_point_masses,
    index_point_charges,
    index_plasma,
    unified_index_field,
    index_from_density,
    
    # Flow field functions
    flow_rotational,
    flow_translational_mass_current,
    flow_charge_current,
    flow_from_vector_potential,
    unified_flow_field,
    magnetic_vector_potential_wire,
    
    # Ray kinematics
    ray_direction_update,
    ray_advection,
    
    # Analysis functions
    estimate_A_M_from_GR,
    current_dragging_ratio,
    charge_to_mass_lensing_ratio,
    pm_gravitomagnetic_field,

    # Massive particle dynamics
    massive_accel_medium,
    newtonian_accel_sum,
    
    # Legacy GR comparison functions
    pm_deflection_angle_point_mass,
    lense_thirring_precession,
    # GW ringdown / echoes
    pm_photon_sphere_radius,
    pm_surface_echo_delay,
    pm_surface_mode_frequency,
    pm_surface_mode_damping_time,
    # Perihelion precession
    pm_perihelion_precession,
    pm_precession_arcsec_per_century,
    pm_integrate_orbit,
)

__all__ = [
    'G', 'c', 'mu_0', 'epsilon_0', 'e', 'm_e',
    'index_point_masses',
    'index_point_charges',
    'index_plasma',
    'unified_index_field',
    'index_from_density',
    'flow_rotational',
    'flow_translational_mass_current',
    'flow_charge_current',
    'flow_from_vector_potential',
    'unified_flow_field',
    'magnetic_vector_potential_wire',
    'ray_direction_update',
    'ray_advection',
    'estimate_A_M_from_GR',
    'current_dragging_ratio',
    'charge_to_mass_lensing_ratio',
    'pm_gravitomagnetic_field',
    'pm_deflection_angle_point_mass',
    'lense_thirring_precession',
    'massive_accel_medium',
    'newtonian_accel_sum',
    'pm_photon_sphere_radius',
    'pm_surface_echo_delay',
    'pm_surface_mode_frequency',
    'pm_surface_mode_damping_time',
    'pm_perihelion_precession',
    'pm_precession_arcsec_per_century',
    'pm_integrate_orbit',
]
