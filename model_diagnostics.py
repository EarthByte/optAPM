"""
Post-run diagnostics for an optimised absolute plate motion model.

Computes, for every timestep of an optimisation run:

  * Net rotation (NR) of the optimised model: the rate (deg/Myr) at 1 Myr sub-steps
    within each timestep interval, summarised per timestep as the median and the
    median absolute deviation (MAD) of those sub-step rates.

  * Trench-normal migration of all subduction zone segments: per timestep, the mean
    and MAD of the trench-orthogonal velocities (mm/yr, positive = retreat/rollback)
    across all trench segments, plus the percentage of segments retreating.

Outputs (in 'model_output/'):

  * <model_name>_diagnostics.csv                - the statistics through time
  * <model_name>_net_rotation.png               - NR median +/- MAD vs time
  * <model_name>_trench_migration.png           - TM mean +/- MAD vs time

These diagnostics run automatically at the end of 'Optimised_APM.py' (if
'generate_diagnostics' is True in 'Optimised_config.py'), and can also be run
standalone on an existing optimised model:

    python model_diagnostics.py

(using the current settings in 'Optimised_config.py', including any OPTAPM_*
environment variable overrides such as OPTAPM_VARIANT / OPTAPM_START_AGE).
"""

import math
import os.path
import sys

import numpy as np
import pygplates

import optimisation_methods
import subduction_convergence_for_absolute_plate_motion as scap


def mad(values, centre):
    """Median absolute deviation about 'centre'."""
    values = np.asarray(values)
    if values.size == 0:
        return float('nan')
    return float(np.median(np.abs(values - centre)))


def calculate_diagnostics(
        data_dir,
        optimised_rotation_filename,    # Relative to the 'data/' directory.
        no_net_rotation_filename,       # Relative to the 'data/' directory.
        trench_resolver,                # TrenchResolver (used to resolve trenches at each time).
        age_range,                      # Iterable of timestep start ages (eg, [5, 10, ..., 1000]).
        interval,
        reference_params_function):     # Function (accepting age) returning (ref_plate_id, _).
    """
    Calculate net rotation and trench migration statistics at each timestep.
    Returns a list of dict rows (one per timestep age).
    """

    optimised_rotation_model = pygplates.RotationModel(
            os.path.join(data_dir, optimised_rotation_filename))
    nn_rotation_model = pygplates.RotationModel(
            os.path.join(data_dir, no_net_rotation_filename))

    rows = []

    for age in age_range:

        ref_rotation_plate_id, _ = reference_params_function(age)

        #
        # Net rotation: rates at 1 Myr sub-steps within [age - interval, age].
        #
        nr_timesteps = range(age - interval, age, 1)
        _, _, PTangle1, _, _, _, _, _, _ = optimisation_methods.ApproximateNR_from_features(
                optimised_rotation_model, nn_rotation_model, nr_timesteps, ref_rotation_plate_id)
        nr_rates = np.abs(np.asarray(PTangle1, dtype=float))  # deg/Myr at each 1 Myr sub-step
        nr_median = float(np.median(nr_rates))
        nr_mad = mad(nr_rates, nr_median)

        #
        # Trench migration: trench-orthogonal velocities across all subduction segments.
        # (Positive = retreat/rollback toward the ocean basin; negative = advance.)
        #
        trench_resolver.generate_resolved_trenches(age)
        trench_features = pygplates.FeatureCollection(
                os.path.join(data_dir, trench_resolver.get_trench_migration_filename()))

        tm_stats = scap.subduction_absolute_motion(
                optimised_rotation_model, trench_features, np.radians(1.0), age,
                velocity_delta_time=interval)
        trench_vel = np.array([stat[2] for stat in tm_stats]) * 10  # cm/yr -> mm/yr
        trench_obl = np.array([stat[3] for stat in tm_stats])
        tm_vel_orth = np.abs(trench_vel) * -np.cos(np.radians(trench_obl))

        tm_mean = float(np.mean(tm_vel_orth))
        tm_median = float(np.median(tm_vel_orth))
        tm_mad = mad(tm_vel_orth, tm_median)
        pct_retreating = float(100.0 * np.mean(tm_vel_orth > 0))

        rows.append({
            'age': age,
            'net_rotation_median_deg_per_myr': nr_median,
            'net_rotation_mad_deg_per_myr': nr_mad,
            'trench_migration_mean_mm_per_yr': tm_mean,
            'trench_migration_median_mm_per_yr': tm_median,
            'trench_migration_mad_mm_per_yr': tm_mad,
            'percent_trenches_retreating': pct_retreating,
        })

        print('  diagnostics at {0} Ma: NR median {1:.3f} deg/Myr, TM mean {2:+.1f} mm/yr ({3:.0f}% retreating)'.format(
                age, nr_median, tm_mean, pct_retreating))
        sys.stdout.flush()

    return rows


def write_diagnostics_csv(rows, csv_filename):
    fields = ['age',
              'net_rotation_median_deg_per_myr', 'net_rotation_mad_deg_per_myr',
              'trench_migration_mean_mm_per_yr', 'trench_migration_median_mm_per_yr',
              'trench_migration_mad_mm_per_yr', 'percent_trenches_retreating']
    with open(csv_filename, 'w') as f:
        f.write(','.join(fields) + '\n')
        for row in rows:
            f.write(','.join('{0:.4f}'.format(row[field]) if field != 'age' else str(row[field])
                             for field in fields) + '\n')
    print('Wrote {0}'.format(csv_filename))


def plot_diagnostics(rows, model_name, output_dir):
    """Generate the net rotation and trench migration plots (PNG)."""
    # Use a non-interactive backend so this works on HPC nodes without a display.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ages = np.array([row['age'] for row in rows])

    #
    # Net rotation plot.
    #
    nr_median = np.array([row['net_rotation_median_deg_per_myr'] for row in rows])
    nr_mad = np.array([row['net_rotation_mad_deg_per_myr'] for row in rows])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.fill_between(ages, nr_median - nr_mad, nr_median + nr_mad,
                    alpha=0.3, label='median +/- MAD (within timestep)')
    ax.plot(ages, nr_median, lw=1.5, marker='o', ms=2.5, label='median')
    # Reference lines: preferred geodynamic upper limit and zero.
    ax.axhline(0.20, ls='--', lw=0.8, color='grey',
               label='0.20 deg/Myr (preferred upper bound)')
    ax.axhline(0.26, ls=':', lw=0.8, color='grey',
               label='0.26 deg/Myr (Conrad & Behn 2010 limit)')
    whole_median = float(np.median(nr_median))
    whole_mad = mad(nr_median, whole_median)
    ax.set_xlabel('Age (Ma)')
    ax.set_ylabel('Net rotation (deg/Myr)')
    ax.set_title('{0}: net rotation through time\n'
                 '(whole-run median {1:.3f} +/- {2:.3f} deg/Myr MAD)'.format(
                    model_name, whole_median, whole_mad), fontsize=10)
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    nr_png = os.path.join(output_dir, '{0}_net_rotation.png'.format(model_name))
    fig.savefig(nr_png, dpi=150)
    plt.close(fig)
    print('Wrote {0}'.format(nr_png))

    #
    # Trench migration plot.
    #
    tm_mean = np.array([row['trench_migration_mean_mm_per_yr'] for row in rows])
    tm_mad = np.array([row['trench_migration_mad_mm_per_yr'] for row in rows])
    pct_retreat = np.array([row['percent_trenches_retreating'] for row in rows])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.fill_between(ages, tm_mean - tm_mad, tm_mean + tm_mad,
                    alpha=0.3, label='mean +/- MAD (across trench segments)')
    ax.plot(ages, tm_mean, lw=1.5, marker='o', ms=2.5, label='mean')
    ax.axhline(0.0, ls='-', lw=0.8, color='black')
    # Observed present-day range of mean trench-normal retreat (Schellart et al. 2008).
    ax.axhspan(9, 15, alpha=0.15, color='green',
               label='observed present-day mean retreat\n(+0.9 to +1.5 cm/yr, Schellart et al. 2008)')
    whole_mean = float(np.mean(tm_mean))
    ax.set_xlabel('Age (Ma)')
    ax.set_ylabel('Trench-normal migration (mm/yr; positive = retreat)')
    ax.set_title('{0}: subduction zone migration through time\n'
                 '(whole-run mean {1:+.1f} mm/yr; mean {2:.0f}% of segments retreating)'.format(
                    model_name, whole_mean, float(np.mean(pct_retreat))), fontsize=10)
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    tm_png = os.path.join(output_dir, '{0}_trench_migration.png'.format(model_name))
    fig.savefig(tm_png, dpi=150)
    plt.close(fig)
    print('Wrote {0}'.format(tm_png))


def generate_model_diagnostics(
        data_dir,
        optimised_rotation_filename,
        no_net_rotation_filename,
        trench_resolver,
        age_range,
        interval,
        reference_params_function,
        model_name,
        output_dir='model_output',
        also_diagnose_nnr_model=True):
    """
    Compute the diagnostics, write the CSV and generate the two plots.
    This is called at the end of a standard optimisation run (see 'Optimised_APM.py').

    If 'also_diagnose_nnr_model' is True then the same diagnostics are also generated
    for the no-net-rotation model itself (with outputs named '<model_name>_NNR_*').
    This shows how subduction zone migration behaves in an NNR world - the zero-NR
    end-member of the uncertainty envelope - and the NNR net rotation plot doubles
    as a sanity check (it should be approximately zero at all times).
    """
    print('Generating model diagnostics (net rotation and trench migration through time)...')
    sys.stdout.flush()

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    rows = calculate_diagnostics(
            data_dir, optimised_rotation_filename, no_net_rotation_filename,
            trench_resolver, age_range, interval, reference_params_function)
    write_diagnostics_csv(rows, os.path.join(output_dir, '{0}_diagnostics.csv'.format(model_name)))
    plot_diagnostics(rows, model_name, output_dir)

    if also_diagnose_nnr_model:
        print('Generating diagnostics for the no-net-rotation model...')
        sys.stdout.flush()
        # Use the NNR model as the rotation model being diagnosed. Its net rotation
        # (relative to itself) should be ~zero - a useful sanity check - while its
        # trench migration statistics show subduction zone kinematics in an NNR world.
        nnr_rows = calculate_diagnostics(
                data_dir, no_net_rotation_filename, no_net_rotation_filename,
                trench_resolver, age_range, interval, reference_params_function)
        nnr_name = '{0}_NNR'.format(model_name)
        write_diagnostics_csv(nnr_rows, os.path.join(output_dir, '{0}_diagnostics.csv'.format(nnr_name)))
        plot_diagnostics(nnr_rows, nnr_name, output_dir)

    return rows


if __name__ == '__main__':

    # Standalone mode: run the diagnostics on an existing optimised model using
    # the current settings in 'Optimised_config.py'.
    from Optimised_config import (datadir, data_model, model_name, start_age, actual_end_age,
                                  interval, get_reference_params, original_rotation_filenames,
                                  topology_filenames)
    from trench_resolver import TrenchResolver

    optimised_rotation_filename = os.path.join(
            data_model, 'optimisation', 'optimised_rotation_model_' + model_name + '.rot')
    no_net_rotation_filename = os.path.join(
            data_model, 'optimisation', 'no_net_rotation_model_' + model_name + '.rot')

    print('Loading topologies...')
    topology_features = []
    for topology_filename in topology_filenames:
        topology_features.extend(pygplates.FeatureCollection(os.path.join(datadir, topology_filename)))

    trench_resolver = TrenchResolver(datadir, original_rotation_filenames, topology_features, data_model)

    age_range = range(actual_end_age + interval, start_age + interval, interval)

    generate_model_diagnostics(
            datadir, optimised_rotation_filename, no_net_rotation_filename,
            trench_resolver, age_range, interval, get_reference_params, model_name)
