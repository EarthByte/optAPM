"""
make_pmag_seed.py - generate a paleomagnetic reference-rotation SEED file for optAPM.

Some plate models store their paleomagnetic ("pmag") absolute path as a separate plate that you
*anchor* to switch the whole model into the pmag reference frame - e.g. the Zahirovic et al. (2022)
model carries the GAPWAP apparent-polar-wander path on plate 701701 (anchor 701701 = PMAG frame).

optAPM seeds each interval's search from `get_rotation(start_age, ref_plate, anchor=0)` of a reference
rotation file (see get_reference_params / OPTAPM_REF_ROTATION_FILE). To seed from such a model's pmag
path you therefore need that path written as a standalone <ref_plate>-relative-to-000 rotation file.

This script extracts the pmag plate's rotation *relative to the reference plate* (the embedded
apparent-polar-wander rotation - identity at present day) and writes it as the reference plate relative
to 000. Then point optAPM at the result:

    OPTAPM_REF_ROTATION_FILE=<path relative to data/> ./run_optapm.sh <model> ...

By default the input rotations, age range and data model are taken from Optimised_config.py, so for the
currently-selected model this is a one-liner:

    python make_pmag_seed.py --pmag-plate 701701

Full options:
    python make_pmag_seed.py --pmag-plate 701701 [--ref-plate 701] [--output <path>]
                             [--data-model NAME] [--rotation-files f1 f2 ...]
                             [--start-age N] [--end-age 0] [--interval N]
"""

import argparse
import math
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pmag-plate', type=int, required=True,
                        help='Plate carrying the paleomagnetic path (the plate you would anchor to get '
                             'the PMAG frame), e.g. 701701 for Zahirovic et al. 2022.')
    parser.add_argument('--ref-plate', type=int, default=701,
                        help='Reference plate the pmag plate is defined relative to, and the plate the '
                             'seed is written for (default 701, Africa).')
    parser.add_argument('--output', default=None,
                        help='Output .rot path (default: data/<data_model>/pmag/<data_model>_pmag_seed.rot).')
    parser.add_argument('--data-model', default=None,
                        help='Data model name (default: the one selected in Optimised_config.py).')
    parser.add_argument('--rotation-files', nargs='+', default=None,
                        help='Input rotation files (default: the input rotations of the data model).')
    parser.add_argument('--start-age', type=float, default=None, help='Oldest age, Ma (default: config).')
    parser.add_argument('--end-age', type=float, default=0.0, help='Youngest age, Ma (default 0).')
    parser.add_argument('--interval', type=float, default=None, help='Sampling step, Ma (default: config).')
    args = parser.parse_args()

    # Let --data-model influence the config defaults (config reads OPTAPM_DATA_MODEL at import).
    if args.data_model:
        os.environ['OPTAPM_DATA_MODEL'] = args.data_model

    import pygplates
    import Optimised_config as cfg

    data_model = args.data_model or cfg.data_model
    datadir = cfg.datadir
    start_age = args.start_age if args.start_age is not None else float(cfg.start_age)
    interval = args.interval if args.interval is not None else float(cfg.interval)
    end_age = args.end_age

    # Input rotation files (absolute paths).
    if args.rotation_files:
        rotation_files = [f if os.path.isabs(f) else os.path.join(datadir, f) for f in args.rotation_files]
    else:
        rotation_files = [os.path.join(datadir, f) for f in cfg.original_rotation_filenames]
    rotation_model = pygplates.RotationModel(rotation_files)

    # Output path.
    if args.output:
        output_path = args.output
    else:
        output_path = os.path.join(datadir, data_model, 'pmag', '{0}_pmag_seed.rot'.format(data_model))
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # Sample R(pmag_plate relative to ref_plate) over [end_age, start_age], written as ref_plate rel 000.
    comment = ('!pmag seed: plate {0} relative to {1} (apparent polar wander) from {2}'
               .format(args.pmag_plate, args.ref_plate, data_model))
    n = int(round((start_age - end_age) / interval))
    times = [end_age + i * interval for i in range(n + 1)]
    lines = []
    for t in times:
        rot = rotation_model.get_rotation(float(t), args.pmag_plate, fixed_plate_id=args.ref_plate)
        if rot.represents_identity_rotation():
            lat, lon, angle = 90.0, 0.0, 0.0
        else:
            lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
        lines.append('{0:<5d} {1:8.1f} {2:10.4f} {3:10.4f} {4:9.4f}  000 {5}\n'.format(
            args.ref_plate, float(t), lat, lon, angle, comment))
    with open(output_path, 'w') as f:
        f.writelines(lines)

    rel_to_data = os.path.relpath(output_path, datadir)
    print('Wrote {0} ({1} samples, {2}-{3} Ma).'.format(output_path, len(lines), int(end_age), int(start_age)))
    print('Use it as the optAPM search seed with:')
    print('    OPTAPM_REF_PLATE_ID={0} OPTAPM_REF_ROTATION_FILE={1} ./run_optapm.sh {2} <cores>'.format(
        args.ref_plate, rel_to_data, data_model))


if __name__ == '__main__':
    # Run from the directory containing Optimised_config.py.
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    main()
