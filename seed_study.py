"""
Seed-count convergence study and profiling harness for the optAPM workflow.

The study answers: how many start seeds (multi-start optimisations) per timestep
are actually needed to reliably find the global minimum of the objective function,
and how much cheaper is a "screen-then-polish" strategy (evaluate the objective
once at every seed, then fully optimise only the best few)?

The harness is checkpointed/resumable so it can also run in environments with
short execution windows: every command can be re-run repeatedly and will continue
where it left off.

Usage (from the optAPM directory; config comes from Optimised_config.py):

    # 1. Prepare data for a timestep (re-run until it prints 'PREP DONE'):
    python seed_study.py prep --age 100

    # 2. Time a single objective evaluation and a single full optimisation:
    python seed_study.py probe --age 100

    # 3. Run the study tasks (re-run until it prints 'ALL TASKS DONE'):
    python seed_study.py run --age 100 [--budget-seconds 1e9]

    # 4. Summarise results:
    python seed_study.py report --age 100

Outputs CSV summary to model_output/seed_study/<age>Ma_report.csv.
"""

import argparse
import json
import marshal
import math
import multiprocessing
import os
import os.path
import pickle
import sys
import time
import types

import numpy as np
import pygplates as pgp

import pmagpy.pmag as pmag

import Optimised_config as cfg


STATE_DIR_TEMPLATE = os.path.join('model_output', 'seed_study')

# Seed-grid sizes to test (requested number of models; for a 'Uniform' search the
# actual number is rounded up to the nearest square).
DEFAULT_SEED_COUNTS = [16, 25, 49, 100, 200]

# Numbers of best screened seeds to polish (for screen-then-polish strategy).
DEFAULT_TOP_NS = [4, 8, 16, 32]


def state_paths(age):
    state_dir = STATE_DIR_TEMPLATE
    if not os.path.isdir(state_dir):
        os.makedirs(state_dir)
    return (os.path.join(state_dir, '{0}Ma_progress.json'.format(age)),
            os.path.join(state_dir, '{0}Ma_state.pkl'.format(age)),
            os.path.join(state_dir, '{0}Ma_tasks.pkl'.format(age)),
            os.path.join(state_dir, '{0}Ma_report.csv'.format(age)))


def load_json(path, default):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def save_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f)
    os.replace(tmp, path)


def load_pickle(path, default=None):
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    return default


def save_pickle(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# Preparation (chunked, resumable)
# ---------------------------------------------------------------------------

def cmd_prep(age, budget_seconds):
    from optapm import ModelSetup as ms
    from no_net_rotation_model import NoNetRotationModel
    from optimised_rotation_updater import OptimisedRotationUpdater
    from plate_velocity_partitioner import PlateVelocityPartitioner
    from continent_fragmentation import ContinentFragmentation
    from trench_resolver import TrenchResolver

    t_start = time.time()
    progress_path, state_path, _, _ = state_paths(age)
    progress = load_json(progress_path, {'nnr_age': None, 'rotfile': False, 'files': False, 'state': False})

    interval = cfg.interval
    ref_rotation_start_age = age
    ref_rotation_end_age = age - interval
    age_range = range(cfg.actual_end_age + interval, cfg.start_age + interval, interval)

    optimisation_sub_dir = os.path.join(cfg.datadir, cfg.data_model, 'optimisation')
    if not os.path.isdir(optimisation_sub_dir):
        os.makedirs(optimisation_sub_dir)

    print('Loading topologies...')
    sys.stdout.flush()
    topology_features = []
    for topology_filename in cfg.topology_filenames:
        topology_features.extend(pgp.FeatureCollection(os.path.join(cfg.datadir, topology_filename)))
    print('  ...{0:.1f} s'.format(time.time() - t_start))
    sys.stdout.flush()

    # --- Optimised rotation file (created once; re-used on later calls) ----
    # Pass end_age != actual_end_age to re-use existing file (continuation mode).
    rot_end_age = age if progress['rotfile'] else cfg.actual_end_age
    optimised_rotation_updater = OptimisedRotationUpdater(
            cfg.datadir, cfg.original_rotation_filenames,
            cfg.start_age, rot_end_age, cfg.actual_end_age, interval,
            cfg.get_reference_params, cfg.data_model, cfg.model_name)
    rotfile = optimised_rotation_updater.get_optimised_rotation_filename()
    if not progress['rotfile']:
        progress['rotfile'] = True
        save_json(progress_path, progress)
    print('Optimised rotation file ready ({0:.1f} s elapsed)'.format(time.time() - t_start))
    sys.stdout.flush()

    # --- No-net-rotation model (chunked) ------------------------------------
    nnr_done_age = progress['nnr_age'] or cfg.actual_end_age
    no_net_rotation_model = NoNetRotationModel(
            cfg.datadir, cfg.original_rotation_filenames, topology_features,
            cfg.start_age, nnr_done_age, cfg.actual_end_age,
            cfg.data_model, cfg.model_name)
    nnr_rotfile = no_net_rotation_model.get_no_net_rotation_filename()

    while nnr_done_age < age:
        if time.time() - t_start > budget_seconds:
            progress['nnr_age'] = nnr_done_age
            save_json(progress_path, progress)
            print('PREP INCOMPLETE (no-net-rotation model at {0}/{1} Ma) - re-run to continue'.format(nnr_done_age, age))
            return
        next_age = min(nnr_done_age + 5, age)
        no_net_rotation_model.update_no_net_rotation(next_age)
        nnr_done_age = next_age
        print('  no-net-rotation model updated to {0} Ma ({1:.1f} s elapsed)'.format(nnr_done_age, time.time() - t_start))
        sys.stdout.flush()
    progress['nnr_age'] = nnr_done_age
    save_json(progress_path, progress)

    # --- Trench + plate velocity files at this age ---------------------------
    print('Resolving trenches and plate velocity points...')
    sys.stdout.flush()
    trench_resolver = TrenchResolver(
            cfg.datadir, cfg.original_rotation_filenames, topology_features, cfg.data_model)
    tm_file = trench_resolver.get_trench_migration_filename()
    trench_resolver.generate_resolved_trenches(ref_rotation_start_age)

    if cfg.plate_velocity_continental_polygons_file:
        plate_velocity_plate_features = list(
            pgp.FeatureCollection(os.path.join(cfg.datadir, cfg.plate_velocity_continental_polygons_file)))
        plate_velocity_features_are_topologies = False
        plate_velocity_fragmentation = ContinentFragmentation(
                cfg.datadir, cfg.original_rotation_filenames, plate_velocity_plate_features,
                cfg.plate_velocity_continental_fragmentation_point_spacing_degrees,
                cfg.plate_velocity_continental_fragmentation_area_threshold_steradians,
                cfg.plate_velocity_continental_fragmentation_gap_threshold_radians,
                age_range)
    else:
        plate_velocity_plate_features = topology_features
        plate_velocity_features_are_topologies = True
        plate_velocity_fragmentation = None

    plate_velocity_partitioner = PlateVelocityPartitioner(
            cfg.datadir, cfg.original_rotation_filenames, plate_velocity_plate_features,
            plate_velocity_features_are_topologies, plate_velocity_fragmentation,
            cfg.data_model, cfg.plate_velocity_grid_spacing)
    pv_file = plate_velocity_partitioner.get_plate_velocity_filename()
    plate_velocity_partitioner.generate_points_and_plate_ids(ref_rotation_start_age)
    print('  ...done ({0:.1f} s elapsed)'.format(time.time() - t_start))
    sys.stdout.flush()

    # --- Reference params, seed current interval of optimised rotation ------
    ref_rotation_plate_id, ref_rotation_file = cfg.get_reference_params(ref_rotation_start_age)
    if ref_rotation_file == cfg.USE_NNR_REFERENCE_FRAME:
        ref_rotation_file = nnr_rotfile
    elif ref_rotation_file == cfg.USE_OPTIMISED_REFERENCE_FRAME:
        ref_rotation_file = rotfile

    _rotation_model = pgp.RotationModel(os.path.join(cfg.datadir, rotfile))
    plate_rotation_005_rel_000 = _rotation_model.get_rotation(ref_rotation_end_age, 5, fixed_plate_id=0)
    plate_rotation_ref_plate_rel_005 = _rotation_model.get_rotation(
            ref_rotation_start_age, ref_rotation_plate_id, fixed_plate_id=5)
    optimised_rotation_updater.update_optimised_rotation(
            plate_rotation_005_rel_000 * plate_rotation_ref_plate_rel_005,
            ref_rotation_plate_id, ref_rotation_start_age)

    # --- Component params, data load, starting conditions --------------------
    enable_fz, fz_w, fz_cf, fz_b = cfg.get_fracture_zone_params(ref_rotation_start_age)
    enable_nr, nr_w, nr_cf, nr_b = cfg.get_net_rotation_params(ref_rotation_start_age)
    enable_tm, tm_w, tm_cf, tm_b = cfg.get_trench_migration_params(ref_rotation_start_age)
    enable_hs, hs_w, hs_cf, hs_b = cfg.get_hotspot_trail_params(ref_rotation_start_age)
    enable_pv, pv_w, pv_cf, pv_b = cfg.get_plate_velocity_params(ref_rotation_start_age)

    params = [cfg.search_radius, cfg.rotation_uncertainty, cfg.search_type, cfg.models,
              cfg.model_stop_condition, cfg.max_iter,
              ref_rotation_plate_id, ref_rotation_start_age, ref_rotation_end_age,
              cfg.interpolation_resolution, cfg.rotation_age_of_interest,
              enable_fz, enable_nr, enable_tm, enable_hs, enable_pv,
              fz_w, nr_w, tm_w, hs_w, pv_w,
              fz_cf, nr_cf, tm_cf, hs_cf, pv_cf,
              fz_b, nr_b, tm_b, hs_b, pv_b,
              cfg.no_auto_ref_rot_longitude, cfg.no_auto_ref_rot_latitude, cfg.no_auto_ref_rot_angle,
              cfg.auto_calc_ref_pole, cfg.search,
              cfg.include_chains, cfg.interpolated_hotspot_trails, cfg.tm_method]

    print('Loading data and computing starting conditions...')
    sys.stdout.flush()
    data = ms.dataLoader(cfg.datadir, rotfile, ref_rotation_file,
                         tm_file=tm_file, pv_file=pv_file, nnr_rotfile=nnr_rotfile,
                         ridge_file=cfg.ridge_file if enable_fz else None,
                         isochron_file=cfg.isochron_file if enable_fz else None,
                         isocob_file=cfg.isocob_file if enable_fz else None,
                         hst_file=cfg.hst_file, hs_file=cfg.hs_file,
                         interpolated_hotspots=cfg.interpolated_hotspots)
    sc = ms.modelStartConditions(params, data, False)

    (x, opt_n, N, lb, ub, model_stop_condition, max_iter,
     rotation_file, _, _, _,
     Lats, Lons, spreading_directions, spreading_asymmetries, seafloor_ages, PID, CPID,
     data_array, weights_array, cost_func_array, bounds_array,
     trench_migration_file, plate_velocity_file, no_net_rotation_file, reformArray, trail_data,
     start_seeds, rotation_age_of_interest_age, data_array_labels_short,
     ref_rot_longitude, ref_rot_latitude, ref_rot_angle, seed_lons, seed_lats) = sc[:35]

    state = {
        'interval': cfg.interval,
        'rotation_file': rotation_file,
        'no_net_rotation_file': no_net_rotation_file,
        'ref_rotation_start_age': ref_rotation_start_age,
        'ref_rotation_end_age': ref_rotation_end_age,
        'ref_rotation_plate_id': ref_rotation_plate_id,
        'Lats': Lats, 'Lons': Lons,
        'spreading_directions': spreading_directions,
        'spreading_asymmetries': spreading_asymmetries,
        'seafloor_ages': seafloor_ages, 'PID': PID, 'CPID': CPID,
        'data_array': data_array, 'weights_array': weights_array,
        'cost_func_code_strings': [marshal.dumps(cf.__code__) for cf in cost_func_array],
        'bounds_array': bounds_array,
        'trench_migration_file': trench_migration_file,
        'plate_velocity_file': plate_velocity_file,
        'reformArray': reformArray, 'trail_data': trail_data,
        'use_trail_age_uncertainty': cfg.use_trail_age_uncertainty,
        'trail_age_uncertainty_ellipse': cfg.trail_age_uncertainty_ellipse,
        'tm_method': cfg.tm_method,
        'x': [list(map(float, xi)) for xi in x],
        'opt_n': opt_n, 'lb': list(lb), 'ub': list(ub),
        'model_stop_condition': model_stop_condition, 'max_iter': max_iter,
        'ref_rot_longitude': float(ref_rot_longitude),
        'ref_rot_latitude': float(ref_rot_latitude),
        'ref_rot_angle': float(ref_rot_angle),
        'search_radius': cfg.search_radius,
    }
    save_pickle(state_path, state)
    progress['files'] = True
    progress['state'] = True
    save_json(progress_path, progress)
    print('PREP DONE ({0:.1f} s total, {1} seeds)'.format(time.time() - t_start, len(x)))


# ---------------------------------------------------------------------------
# Worker processes
# ---------------------------------------------------------------------------

_g_obj_f = None
_g_nlopt_args = None


def _init_worker(state):
    global _g_obj_f, _g_nlopt_args
    from objective_function import ObjectiveFunction
    cost_func_array = [types.FunctionType(marshal.loads(code), globals(), 'cost_func')
                       for code in state['cost_func_code_strings']]
    _g_obj_f = ObjectiveFunction(
            state['interval'], state['rotation_file'], state['no_net_rotation_file'],
            state['ref_rotation_start_age'], state['Lats'], state['Lons'],
            state['spreading_directions'], state['spreading_asymmetries'], state['seafloor_ages'],
            state['PID'], state['CPID'], state['data_array'], state['weights_array'],
            cost_func_array, state['bounds_array'],
            state['trench_migration_file'], state['plate_velocity_file'],
            state['ref_rotation_end_age'], state['ref_rotation_plate_id'],
            state['reformArray'], state['trail_data'],
            state['use_trail_age_uncertainty'], state['trail_age_uncertainty_ellipse'],
            state['tm_method'])
    _g_nlopt_args = (state['opt_n'], state['lb'], state['ub'],
                     state['model_stop_condition'], state['max_iter'])


def _run_task(task):
    kind = task['kind']
    seed = task['seed']
    t0 = time.time()
    if kind == 'screen':
        cost = _g_obj_f(seed, None)
        return task['id'], {'cost': float(cost), 'time': time.time() - t0}
    else:  # full optimisation
        import nlopt
        opt_n, lb, ub, model_stop_condition, max_iter = _g_nlopt_args
        opt = nlopt.opt(nlopt.LN_COBYLA, opt_n)
        opt.set_min_objective(_g_obj_f)
        opt.set_lower_bounds(lb)
        opt.set_upper_bounds(ub)
        if model_stop_condition != 'threshold':
            opt.set_maxeval(max_iter)
        else:
            opt.set_ftol_rel(1e-6)
            opt.set_xtol_rel(1e-8)
            # Safety cap: COBYLA occasionally fails to converge (cycles for thousands of
            # evaluations). The vast majority of optimisations converge in < 200 evaluations,
            # so 300 is a generous cap that only affects pathological cases.
            opt.set_maxeval(300)
        xopt = opt.optimize(list(seed))
        return task['id'], {'xopt': [float(v) for v in xopt],
                            'cost': float(opt.last_optimum_value()),
                            'nevals': opt.get_numevals(),
                            'time': time.time() - t0}


# ---------------------------------------------------------------------------
# Probe
# ---------------------------------------------------------------------------

def cmd_probe(age):
    _, state_path, _, _ = state_paths(age)
    state = load_pickle(state_path)
    if state is None:
        print('No state found - run prep first.')
        return
    _init_worker(state)
    x = state['x']
    # Time single evaluations.
    times = []
    for seed in x[:3]:
        t0 = time.time()
        cost = _g_obj_f(seed, None)
        times.append(time.time() - t0)
        print('eval at ({0:.1f}, {1:.1f}, {2:.2f}) -> cost {3:.4f} in {4:.3f} s'.format(
                seed[0], seed[1], seed[2], cost, times[-1]))
    # Time one full optimisation.
    _, res = _run_task({'kind': 'full', 'seed': x[0], 'id': 'probe'})
    print('full COBYLA: cost {0:.4f}, {1} evals, {2:.1f} s'.format(res['cost'], res['nevals'], res['time']))
    print('mean eval time: {0:.3f} s'.format(np.mean(times)))


# ---------------------------------------------------------------------------
# Run (chunked, resumable)
# ---------------------------------------------------------------------------

def generate_uniform_seeds(num_models, search_radius, ref_rot_longitude, ref_rot_latitude, ref_rot_angle):
    """Uniform seed grid exactly as ModelSetup.modelStartConditions (search_type='Uniform')."""
    num_points_float = num_models * 1.0 / (0.5 * (1 - math.cos(math.radians(search_radius))))
    n = int(math.ceil(math.sqrt(num_points_float)))
    seeds = []
    for lon_index in range(n):
        for lat_index in range(n):
            theta = 2 * np.pi * ((lon_index + 0.5) / float(n))
            phi = np.arccos(2 * ((lat_index + 0.5) / float(n)) - 1.0)
            xyz = (np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi))
            point = pgp.convert_point_on_sphere_to_lat_lon_point(xyz)
            lat, lon = point.get_latitude(), point.get_longitude()
            if lat > 90 - search_radius:
                rlon, rlat = pmag.dodirot(lon, lat, ref_rot_longitude, ref_rot_latitude)
                seeds.append([float(rlon), float(rlat), float(ref_rot_angle)])
    return seeds


def build_tasks(state, seed_counts, top_ns):
    """Build the full task list (screen + full multistart + reduced grids)."""
    tasks = {}
    x = state['x']
    for i, seed in enumerate(x):
        tasks['screen_{0}'.format(i)] = {'kind': 'screen', 'seed': seed, 'id': 'screen_{0}'.format(i)}
        tasks['full_{0}'.format(i)] = {'kind': 'full', 'seed': seed, 'id': 'full_{0}'.format(i)}
    for n in seed_counts:
        seeds_n = generate_uniform_seeds(n, state['search_radius'],
                                         state['ref_rot_longitude'], state['ref_rot_latitude'],
                                         state['ref_rot_angle'])
        for i, seed in enumerate(seeds_n):
            tid = 'grid{0}_{1}'.format(n, i)
            tasks[tid] = {'kind': 'full', 'seed': seed, 'id': tid}
    # Note: 'polish' tasks for screen-then-polish are just the 'full_<i>' tasks of the
    # top-N screened seeds, so they don't need separate optimisations - the report
    # simply combines screening costs with the relevant full results.
    return tasks


def cmd_run(age, budget_seconds, processes, seed_counts, top_ns):
    t_start = time.time()
    progress_path, state_path, tasks_path, _ = state_paths(age)
    state = load_pickle(state_path)
    if state is None:
        print('No state found - run prep first.')
        return

    task_store = load_pickle(tasks_path)
    if task_store is None:
        tasks = build_tasks(state, seed_counts, top_ns)
        task_store = {'tasks': tasks, 'done': {}, 'seed_counts': seed_counts, 'top_ns': top_ns}
        save_pickle(tasks_path, task_store)

    tasks = task_store['tasks']
    done = task_store['done']
    pending = [tid for tid in tasks if tid not in done]
    # Run screening tasks first (cheap, and needed for screen-then-polish).
    pending.sort(key=lambda tid: (0 if tasks[tid]['kind'] == 'screen' else 1, tid))
    print('{0} tasks pending, {1} done'.format(len(pending), len(done)))
    if not pending:
        print('ALL TASKS DONE')
        return

    pool = multiprocessing.Pool(processes, initializer=_init_worker, initargs=(state,))
    n_completed = 0
    try:
        result_iter = pool.imap_unordered(_run_task, (tasks[tid] for tid in pending))
        for tid, res in result_iter:
            done[tid] = res
            n_completed += 1
            if n_completed % 20 == 0:
                save_pickle(tasks_path, task_store)
            if time.time() - t_start > budget_seconds:
                break
    finally:
        pool.terminate()
        pool.join()
        save_pickle(tasks_path, task_store)

    remaining = len(tasks) - len(done)
    print('completed {0} tasks this call ({1:.1f} s); {2} remaining'.format(
            n_completed, time.time() - t_start, remaining))
    if remaining == 0:
        print('ALL TASKS DONE')


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def cmd_report(age):
    _, state_path, tasks_path, report_path = state_paths(age)
    state = load_pickle(state_path)
    task_store = load_pickle(tasks_path)
    if state is None or task_store is None:
        print('No state/tasks found - run prep and run first.')
        return
    tasks, done = task_store['tasks'], task_store['done']
    seed_counts, top_ns = task_store['seed_counts'], task_store['top_ns']
    x = state['x']
    n_seeds = len(x)

    def pole_dist_deg(a, b):
        import geoTools
        lat1, lon1 = geoTools.checkLatLon(a[1], a[0])
        lat2, lon2 = geoTools.checkLatLon(b[1], b[0])
        p1 = pgp.PointOnSphere(lat1, lon1)
        p2 = pgp.PointOnSphere(lat2, lon2)
        return math.degrees(pgp.GeometryOnSphere.distance(p1, p2))

    rows = []

    # Ground truth: full multistart.
    full_results = [done.get('full_{0}'.format(i)) for i in range(n_seeds)]
    have_full = [r for r in full_results if r]
    if not have_full:
        print('No full multistart results yet.')
        return
    gt = min(have_full, key=lambda r: r['cost'])
    gt_cost, gt_x = gt['cost'], gt['xopt']
    total_full_time = sum(r['time'] for r in have_full)
    total_full_evals = sum(r['nevals'] for r in have_full)
    mean_opt_time = total_full_time / len(have_full)
    print('Ground truth ({0}/{1} full optimisations): cost {2:.6f} at '
          'lon {3:.2f} lat {4:.2f} ang {5:.3f}'.format(
            len(have_full), n_seeds, gt_cost, gt_x[0], gt_x[1], gt_x[2]))
    print('  mean {0:.0f} evals and {1:.1f} s per optimisation; total {2:.0f} s CPU'.format(
            total_full_evals / len(have_full), mean_opt_time, total_full_time))
    rows.append(['full_multistart', n_seeds, len(have_full), '{0:.6f}'.format(gt_cost),
                 '{0:.3f}'.format(gt_x[0]), '{0:.3f}'.format(gt_x[1]), '{0:.4f}'.format(gt_x[2]),
                 '0.0', '0.0', '{0:.1f}'.format(total_full_time)])

    # Reduced grids.
    for n in seed_counts:
        grid_results = [done[tid] for tid in done if tid.startswith('grid{0}_'.format(n))]
        n_grid_tasks = sum(1 for tid in tasks if tid.startswith('grid{0}_'.format(n)))
        if not grid_results or len(grid_results) < n_grid_tasks:
            print('grid {0}: incomplete ({1}/{2})'.format(n, len(grid_results), n_grid_tasks))
            continue
        best = min(grid_results, key=lambda r: r['cost'])
        t_total = sum(r['time'] for r in grid_results)
        cost_diff_pct = 100.0 * (best['cost'] - gt_cost) / abs(gt_cost)
        pdist = pole_dist_deg(best['xopt'], gt_x)
        print('grid {0:4d} ({1:3d} seeds): best cost {2:.6f} (+{3:.3f}%), pole dist {4:6.2f} deg, '
              'ang diff {5:+.3f}, CPU {6:.0f} s'.format(
                n, len(grid_results), best['cost'], cost_diff_pct, pdist,
                best['xopt'][2] - gt_x[2], t_total))
        rows.append(['uniform_grid_{0}'.format(n), n, len(grid_results), '{0:.6f}'.format(best['cost']),
                     '{0:.3f}'.format(best['xopt'][0]), '{0:.3f}'.format(best['xopt'][1]),
                     '{0:.4f}'.format(best['xopt'][2]),
                     '{0:.4f}'.format(cost_diff_pct), '{0:.2f}'.format(pdist), '{0:.1f}'.format(t_total)])

    # Screen-then-polish.
    screen_results = [done.get('screen_{0}'.format(i)) for i in range(n_seeds)]
    if all(screen_results):
        screen_costs = np.array([r['cost'] for r in screen_results])
        t_screen = sum(r['time'] for r in screen_results)
        order = np.argsort(screen_costs)
        for top_n in top_ns:
            top_idx = order[:top_n]
            polish = [done.get('full_{0}'.format(i)) for i in top_idx]
            if not all(polish):
                print('screen+polish top {0}: polish results incomplete'.format(top_n))
                continue
            best = min(polish, key=lambda r: r['cost'])
            t_total = t_screen + sum(r['time'] for r in polish)
            cost_diff_pct = 100.0 * (best['cost'] - gt_cost) / abs(gt_cost)
            pdist = pole_dist_deg(best['xopt'], gt_x)
            print('screen all {0} seeds + polish top {1:3d}: best cost {2:.6f} (+{3:.3f}%), '
                  'pole dist {4:6.2f} deg, CPU {5:.0f} s ({6:.1f}x speedup vs full)'.format(
                    n_seeds, top_n, best['cost'], cost_diff_pct, pdist, t_total,
                    total_full_time / t_total))
            rows.append(['screen_polish_top{0}'.format(top_n), n_seeds, top_n,
                         '{0:.6f}'.format(best['cost']),
                         '{0:.3f}'.format(best['xopt'][0]), '{0:.3f}'.format(best['xopt'][1]),
                         '{0:.4f}'.format(best['xopt'][2]),
                         '{0:.4f}'.format(cost_diff_pct), '{0:.2f}'.format(pdist), '{0:.1f}'.format(t_total)])
    else:
        print('screening incomplete: {0}/{1}'.format(sum(1 for r in screen_results if r), n_seeds))

    with open(report_path, 'w') as f:
        f.write('strategy,seeds_requested,optimisations_run,best_cost,lon,lat,angle,cost_diff_pct,pole_dist_deg,cpu_time_s\n')
        for row in rows:
            f.write(','.join(str(v) for v in row) + '\n')
    print('Wrote {0}'.format(report_path))


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='optAPM seed-count study (resumable)')
    parser.add_argument('command', choices=['prep', 'probe', 'run', 'report'])
    parser.add_argument('--age', type=int, required=True)
    parser.add_argument('--budget-seconds', type=float, default=1e9,
                        help='Approximate wall-time budget for this call (prep/run are resumable)')
    parser.add_argument('--processes', type=int, default=multiprocessing.cpu_count())
    parser.add_argument('--seed-counts', type=int, nargs='+', default=DEFAULT_SEED_COUNTS)
    parser.add_argument('--top-ns', type=int, nargs='+', default=DEFAULT_TOP_NS)
    args = parser.parse_args()

    if args.age % cfg.interval != 0:
        raise SystemExit('age must be a multiple of the interval ({0})'.format(cfg.interval))

    if args.command == 'prep':
        cmd_prep(args.age, args.budget_seconds)
    elif args.command == 'probe':
        cmd_probe(args.age)
    elif args.command == 'run':
        cmd_run(args.age, args.budget_seconds, args.processes, args.seed_counts, args.top_ns)
    elif args.command == 'report':
        cmd_report(args.age)


if __name__ == '__main__':
    main()
