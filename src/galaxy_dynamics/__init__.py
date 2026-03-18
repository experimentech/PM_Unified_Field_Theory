"""Galaxy dynamics package.

Exports rotation curve modeling utilities and mock SPARC data loader.
"""

from .rotation import (
    DiskParams, MediumParams,
    mass_enclosed_exponential, accel_baryonic, delta_n_medium, accel_medium,
    circular_velocity, rotation_curve, deflection_angle_axisymmetric
)
from .data import RotationCurve, load_sparc_mock, load_sparc_real
from .fitting import fit_rotation_curve, chi_square, fit_population, compute_residual_metrics
from .halos import (
    NFWParams, BurkertParams, mass_enclosed_nfw, mass_enclosed_burkert,
    circular_velocity_nfw, circular_velocity_burkert, halo_velocity_profile,
    fit_halo_rotation_curve, fit_disk_halo_rotation_curve
)
from .compare import compare_models, aggregate_statistics, export_comparison_results
from .relations import (
    extract_btf_points, fit_btf, compute_btf,
    extract_rar_points, compute_rar
)

__all__ = [
    'DiskParams', 'MediumParams', 'mass_enclosed_exponential', 'accel_baryonic',
    'delta_n_medium', 'accel_medium', 'circular_velocity', 'rotation_curve',
    'deflection_angle_axisymmetric', 'RotationCurve', 'load_sparc_mock', 'load_sparc_real',
    'fit_rotation_curve', 'chi_square', 'fit_population', 'compute_residual_metrics',
    'NFWParams', 'BurkertParams', 'mass_enclosed_nfw', 'mass_enclosed_burkert',
    'circular_velocity_nfw', 'circular_velocity_burkert', 'halo_velocity_profile',
    'fit_halo_rotation_curve', 'fit_disk_halo_rotation_curve',
    'compare_models', 'aggregate_statistics', 'export_comparison_results',
    'extract_btf_points','fit_btf','compute_btf','extract_rar_points','compute_rar'
]
