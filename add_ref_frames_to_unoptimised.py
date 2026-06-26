"""
add_ref_frames_to_unoptimised.py
================================

Post-processing step of the optAPM workflow.

After the optimisation has been run, this script reads the original un-optimised rotation
files together with the 005-000 reference-frame rotations produced by the workflow, and
writes a new set of rotation files that contain ALL the available absolute reference frames
of a plate model as distinct (high-numbered) plate IDs:

    14  pmag               - the model's native paleomagnetic reference frame
    15  optimised          - the 'best' optimised mantle reference frame
    16  optimised_max_NNR  - the 'nr_max' optimised mantle reference frame
                             (the maximum net-rotation end-member of the uncertainty envelope)
    17  no_net_rotation    - the no-net-rotation reference frame (the zero net-rotation end-member)
    18  true_polar_wander  - an approximation of true polar wander, TPW = pmag - optimised(best)

This is the step that "assembles a comprehensive rotation file with all the available
reference frames for a given plate model". Each optimisation run produces its own optimised
005-000 reference frame, so this script requires both optimisation runs (the default 'best'
run and the 'nr_max' run) plus the no-net-rotation file that every run also produces:

    1) python Optimised_APM.py                       # 'best' run  -> optimised + no-net-rotation
    2) OPTAPM_VARIANT=nr_max python Optimised_APM.py  # 'nr_max' run -> optimised_max_NNR
    3) python add_ref_frames_to_unoptimised.py        # this script -> comprehensive rotation file(s)

The resulting rotation hierarchy looks like:

           015  016  017  018
            \    \   |    /
             \------------/
                   |
                  000          <- the "default" reference frame (one of 014-018)
                   |
                  014          <- pmag (the un-optimised model is attached here)
                   |
               ---------
               |       |
              101     701 ...

Anchoring (in GPlates/pyGPlates) to any of 014/015/016/017/018 reconstructs the model in
that reference frame. Anchoring 000 uses whichever frame `default_reference_frame_plate_id`
is set to below.

IMPORTANT - reference frame of the un-optimised input model. This script requires the
un-optimised input model to be in its paleomagnetic frame (anchoring plate 000 of the input
gives PMAG), because plate 014 is set to that native input frame and the optimised / no-net
rotation 005-000 frames are all expressed relative to it. This requirement is exposed as the
clearly-flagged `unoptimised_input_is_in_pmag_frame` switch in the "Input parameters" section
below; the script aborts if it is not set True. If a model's native 000 frame is NOT
paleomagnetic (for example the Zahirovic et al. 2022 model, whose native frame is a mantle
frame and whose PMAG frame is reached by anchoring plate 701701), the input must first be
re-expressed so that anchor 000 = PMAG before running this script.

By default all input/output filenames are derived from "Optimised_config.py" (the same config
the optimisation used), so this script normally needs no editing. Every parameter can still be
overridden in the "Input parameters" section below.
"""

import os.path
import pygplates

import Optimised_config as config


#### Input parameters ####

# =====================================================================================
# IMPORTANT  -  reference frame of the un-optimised INPUT model  (read this!)
# =====================================================================================
# The un-optimised input rotation files MUST be in their paleomagnetic (PMAG) frame:
# anchoring plate 000 of the input must give PMAG. Plate <pmag_ref_frame_plate_id> (014)
# in the output is set to this native input frame, and ALL of the optimised / no-net
# rotation 005-000 reference frames are expressed relative to it. The whole script (and
# the optimisation that produced the 005-000 files) depends on this assumption.
#
# Set this True ONLY once you have confirmed the input is in its PMAG frame.
#
# If your model's native 000 frame is NOT paleomagnetic - e.g. the Zahirovic et al. (2022)
# model, whose native frame is a mantle frame and whose PMAG frame is reached by anchoring
# plate 701701 - you must FIRST re-express the input rotations so that anchor 000 = PMAG,
# otherwise plate 014 will not be the paleomagnetic frame. (The script aborts below if this
# is left False, rather than silently mislabelling the reference frames.)
# =====================================================================================
unoptimised_input_is_in_pmag_frame = True

# Directory containing the data model (the 'data/<data_model>/' directory). All input and
# output paths below are resolved relative to this directory, so the script can be run from
# anywhere (it does not need to be run from inside the data model directory).
data_model_dir = os.path.join(config.datadir, config.data_model)

# The 'best' (default variant) model name and the 'nr_max' variant model name.
#
# "Optimised_config.py" appends the run variant to the model name for any variant other than
# 'best', so reconstruct the base ('best') name robustly regardless of any OPTAPM_VARIANT in
# the current environment.
_base_model_name = config.model_name
if config.run_variant != 'best' and _base_model_name.endswith('_' + config.run_variant):
    _base_model_name = _base_model_name[:-len('_' + config.run_variant)]
best_model_name = _base_model_name
nr_max_model_name = _base_model_name + '_nr_max'

# Original (un-optimised) rotation files (absolute paths).
#
# Taken from the optimisation config so the same input files are used. These are the files
# whose fixed plate 000 (the PMAG frame) becomes plate 014 in the output.
input_rotation_filenames = [os.path.join(config.datadir, rotation_filename)
                            for rotation_filename in config.original_rotation_filenames]

# The optimised / no-net-rotation 005-000 reference-frame files produced by the two runs.
input_optimised_rotation_file = os.path.join(
    data_model_dir, 'optimisation', 'optimised_rotation_model_' + best_model_name + '.rot')
input_optimised_max_NNR_rotation_file = os.path.join(
    data_model_dir, 'optimisation', 'optimised_rotation_model_' + nr_max_model_name + '.rot')
input_no_net_rotation_file = os.path.join(
    data_model_dir, 'optimisation', 'no_net_rotation_model_' + best_model_name + '.rot')

# The suffix appended to the basename of each input rotation file to form the output filename.
# The output files are written to the 'optimisation/' sub-directory of the data model directory.
output_rotation_filenames_suffix = '_' + best_model_name + '_with_ref_frames'

# Plate IDs of the reference frames.
pmag_ref_frame_plate_id = 14
optimised_ref_frame_plate_id = 15
optimised_max_NNR_ref_frame_plate_id = 16
no_net_rotation_ref_frame_plate_id = 17
true_polar_wander_ref_frame_plate_id = 18

# The reference frame that plate 000 should correspond to (must be one of the IDs above).
default_reference_frame_plate_id = pmag_ref_frame_plate_id

# The plate ID used to express true polar wander (TPW = pmag - optimised). TPW is the residual
# rotation of pmag not captured by the optimised mantle frame; following Tetley et al. (2019)
# it is evaluated on Africa (701).
true_polar_wander_reference_plate_id = 701

##########################


# A short human-readable name for each reference frame plate ID (used in rotation descriptions).
reference_frame_names = {
    pmag_ref_frame_plate_id:              'paleomag',
    optimised_ref_frame_plate_id:         'optimised mantle',
    optimised_max_NNR_ref_frame_plate_id: 'optimised mantle (max net-rotation)',
    no_net_rotation_ref_frame_plate_id:   'no-net rotation',
    true_polar_wander_ref_frame_plate_id: 'approx true polar wander',
}
if default_reference_frame_plate_id not in reference_frame_names:
    raise ValueError('default_reference_frame_plate_id is not one of the reference frame plate IDs')

# Refuse to run unless the user has confirmed the input is in its paleomagnetic frame
# (see the prominent note in the "Input parameters" section above).
if not unoptimised_input_is_in_pmag_frame:
    raise RuntimeError(
        "'unoptimised_input_is_in_pmag_frame' is False: this script requires the un-optimised "
        "input model to be in its paleomagnetic frame (anchoring plate 000 of the input gives "
        "PMAG), because plate {0} is set to that frame and the optimised/NNR frames are expressed "
        "relative to it. Re-express the input so that anchor 000 = PMAG (e.g. the Zahirovic et al. "
        "2022 model's PMAG frame is reached via anchor plate 701701), then set this switch True."
        .format(pmag_ref_frame_plate_id))


#
# Load the original (un-optimised) rotation files and relabel their PMAG frame (fixed plate 000)
# as plate <pmag_ref_frame_plate_id> (014). Each input file gets a corresponding output file.
#
output_rotation_feature_collections = [pygplates.FeatureCollection(filename)
                                       for filename in input_rotation_filenames]

output_rotation_filenames = []
for filename in input_rotation_filenames:
    basename, ext = os.path.splitext(os.path.basename(filename))
    output_rotation_filenames.append(
        os.path.join(data_model_dir, 'optimisation', basename + output_rotation_filenames_suffix + ext))

# Change fixed plate 000 to <pmag_ref_frame_plate_id> (and drop any leftover 005-000 features).
max_time_of_sequences_with_fixed_000 = 0.0
for rotation_feature_collection in output_rotation_feature_collections:

    def is_feature_005_000(feature):
        total_reconstruction_pole = feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, _ = total_reconstruction_pole
            return fixed_plate_id == 0 and moving_plate_id == 5
        return False
    rotation_feature_collection.remove(is_feature_005_000)

    for rotation_feature in rotation_feature_collection:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            if fixed_plate_id == 0:
                # Everything that referenced 000 (the PMAG frame) now references plate 014.
                rotation_feature.set_total_reconstruction_pole(
                    pmag_ref_frame_plate_id, moving_plate_id, rotation_sequence)
                max_time = max(ts.get_time() for ts in rotation_sequence.get_time_samples())
                max_time_of_sequences_with_fixed_000 = max(max_time_of_sequences_with_fixed_000, max_time)

# The un-optimised model anchored to the PMAG frame (plate 014). Used for the TPW calculation.
original_rotation_model_anchored_pmag = pygplates.RotationModel(
    output_rotation_feature_collections,
    default_anchor_plate_id=pmag_ref_frame_plate_id)


#
# Helpers.
#
def find_005_000_rotation_sequence(rotation_filename):
    """Return the (moving 005, fixed 000) rotation sequence stored in an optimised/NNR file."""
    rotation_sequence_005_000 = None
    for rotation_feature in pygplates.FeatureCollection(rotation_filename):
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            if fixed_plate_id == 0 and moving_plate_id == 5:
                if rotation_sequence_005_000 is not None:
                    raise ValueError('Found multiple 005-000 sequences in {0}'.format(rotation_filename))
                rotation_sequence_005_000 = rotation_sequence
    if rotation_sequence_005_000 is None:
        raise ValueError('Found no 005-000 sequence in {0}'.format(rotation_filename))
    return rotation_sequence_005_000


def mantle_frame_samples_relative_to_pmag(rotation_sequence_005_000):
    """
    Convert an optimised/NNR 005-000 sequence into the rotation of that mantle frame relative
    to PMAG, R(pmag -> frame), as a list of (time, finite_rotation) pairs.

    In an optimised output file the 005-000 feature stores R(000->005) where, by construction,
    000 is the mantle frame and 005 is PMAG. Hence R(pmag->frame) = R(005->000) = inverse[R(000->005)].
    """
    return [(ts.get_time(), ts.get_value().get_finite_rotation().get_inverse(), ts.is_enabled())
            for ts in rotation_sequence_005_000]


def true_polar_wander_samples_relative_to_pmag(optimised_rotation_sequence_005_000):
    """
    R(pmag -> TPW) as a list of (time, finite_rotation) pairs.

    TPW is the part of the paleomagnetic frame not captured by the optimised mantle frame
    (TPW = pmag - optimised), evaluated on Africa (701):

        R(pmag->TPW) = R(pmag->701) * R(000->005)_optimised

    where R(pmag->701) is the un-optimised Africa rotation in the PMAG frame and R(000->005)
    is the stored optimised reference-frame rotation.
    """
    samples = []
    for ts in optimised_rotation_sequence_005_000:
        time = ts.get_time()
        optimised_finite_rotation = ts.get_value().get_finite_rotation()           # R(000->005)
        africa_in_pmag = original_rotation_model_anchored_pmag.get_rotation(        # R(pmag->701)
            time, true_polar_wander_reference_plate_id)
        samples.append((time, africa_in_pmag * optimised_finite_rotation, ts.is_enabled()))
    return samples


def identity_samples(begin_time, end_time=0.0):
    """Two identity samples spanning [end_time, begin_time]."""
    return [(end_time, pygplates.FiniteRotation(), True),
            (begin_time, pygplates.FiniteRotation(), True)]


# Read the three optimised/NNR 005-000 reference-frame sequences.
optimised_sequence = find_005_000_rotation_sequence(input_optimised_rotation_file)
optimised_max_NNR_sequence = find_005_000_rotation_sequence(input_optimised_max_NNR_rotation_file)
no_net_rotation_sequence = find_005_000_rotation_sequence(input_no_net_rotation_file)

# Each reference frame expressed RELATIVE TO PMAG (014): R(pmag -> frame) sample lists.
ref_frame_samples_relative_to_pmag = {
    pmag_ref_frame_plate_id:              identity_samples(max_time_of_sequences_with_fixed_000),
    optimised_ref_frame_plate_id:         mantle_frame_samples_relative_to_pmag(optimised_sequence),
    optimised_max_NNR_ref_frame_plate_id: mantle_frame_samples_relative_to_pmag(optimised_max_NNR_sequence),
    no_net_rotation_ref_frame_plate_id:   mantle_frame_samples_relative_to_pmag(no_net_rotation_sequence),
    true_polar_wander_ref_frame_plate_id: true_polar_wander_samples_relative_to_pmag(optimised_sequence),
}

# A rotation model of the reference frames relative to PMAG, used to re-express them relative to
# whichever frame the default (000) is set to: R(default->frame) = inverse[R(pmag->default)] * R(pmag->frame).
ref_frames_relative_to_pmag_features = []
for ref_frame_plate_id, samples in ref_frame_samples_relative_to_pmag.items():
    if ref_frame_plate_id == pmag_ref_frame_plate_id:
        continue  # pmag is the anchor (identity relative to itself); no feature needed.
    ref_frames_relative_to_pmag_features.append(
        pygplates.Feature.create_total_reconstruction_sequence(
            fixed_plate_id=pmag_ref_frame_plate_id,
            moving_plate_id=ref_frame_plate_id,
            total_reconstruction_pole=pygplates.GpmlIrregularSampling(
                [pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(rot), time, '', enabled)
                 for time, rot, enabled in samples])))
ref_frames_relative_to_pmag_model = pygplates.RotationModel(
    ref_frames_relative_to_pmag_features, default_anchor_plate_id=pmag_ref_frame_plate_id)


#
# Build the output reference-frame features, each stored as a (fixed 000, moving <frame>) sequence
# so that anchoring 000 gives the chosen default frame and anchoring <frame> gives that frame.
#
# Stored finite rotation is R(000->frame) = R(default->frame) = inverse[R(pmag->default)] * R(pmag->frame).
#
output_reference_frame_features = []
for ref_frame_plate_id, samples_relative_to_pmag in ref_frame_samples_relative_to_pmag.items():

    name = reference_frame_names[ref_frame_plate_id]

    if ref_frame_plate_id == pmag_ref_frame_plate_id:
        # The PMAG feature carries R(000->014) = R(default->pmag) = inverse[R(pmag->default)].
        # Sample it at the default frame's native times (or identity span if default is PMAG).
        if default_reference_frame_plate_id == pmag_ref_frame_plate_id:
            sample_times = [t for t, _, _ in ref_frame_samples_relative_to_pmag[pmag_ref_frame_plate_id]]
        else:
            sample_times = [t for t, _, _ in ref_frame_samples_relative_to_pmag[default_reference_frame_plate_id]]
        out_samples = []
        for time in sample_times:
            pmag_to_default = ref_frames_relative_to_pmag_model.get_rotation(
                time, default_reference_frame_plate_id)
            out_samples.append((time, pmag_to_default.get_inverse(), True))
    else:
        out_samples = []
        for time, rot_pmag_to_frame, enabled in samples_relative_to_pmag:
            pmag_to_default = ref_frames_relative_to_pmag_model.get_rotation(
                time, default_reference_frame_plate_id)
            out_samples.append((time, pmag_to_default.get_inverse() * rot_pmag_to_frame, enabled))

    if ref_frame_plate_id == default_reference_frame_plate_id:
        description = ' Reference frames: {0} ({1:03d} and 000)'.format(name, ref_frame_plate_id)
    else:
        default_name = reference_frame_names[default_reference_frame_plate_id]
        description = ' Reference frames: {0} ({1:03d}) and {2} (000)'.format(
            name, ref_frame_plate_id, default_name)

    output_reference_frame_features.append(
        pygplates.Feature.create_total_reconstruction_sequence(
            fixed_plate_id=0,
            moving_plate_id=ref_frame_plate_id,
            total_reconstruction_pole=pygplates.GpmlIrregularSampling(
                [pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(rot), time, description, enabled)
                 for time, rot, enabled in out_samples])))


#
# Write the output. The reference-frame features are bundled into the first output rotation file
# (rather than a separate 'reference_frames.rot') so that 000 is always defined wherever the model
# is loaded - a novice user cannot accidentally drop the file that defines plate 000.
#
output_reference_frames_feature_collection = pygplates.FeatureCollection(output_reference_frame_features)
for rotation_feature in output_rotation_feature_collections[0]:
    output_reference_frames_feature_collection.add(rotation_feature)
output_rotation_feature_collections[0] = output_reference_frames_feature_collection

for output_rotation_feature_collection, output_rotation_filename in zip(
        output_rotation_feature_collections, output_rotation_filenames):
    output_rotation_feature_collection.write(output_rotation_filename)
    print('Wrote {0}'.format(output_rotation_filename))

print('Reference frames added: pmag={0}, optimised={1}, optimised_max_NNR={2}, '
      'no_net_rotation={3}, true_polar_wander={4}  (default 000 = {5})'.format(
          pmag_ref_frame_plate_id, optimised_ref_frame_plate_id, optimised_max_NNR_ref_frame_plate_id,
          no_net_rotation_ref_frame_plate_id, true_polar_wander_ref_frame_plate_id,
          reference_frame_names[default_reference_frame_plate_id]))
