from collections import defaultdict
from configparser import ConfigParser
import itertools
from funlib.math import cantor_number
from synister import SynisterDb, find_optimal_split, ImpossibleSplit
import argparse
import json
import numpy as np
import random

DATASETS = {
    "fafb": {
        "files": [
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_v4.csv',
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_verified_predicted_synapses_by_transmitter_v4.csv',
        ],
        "holdout_files": [
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_confident_v4.csv',
        ],
        "voxel_size": (40, 4, 4),
        "db_name": 'synister_fafb_v4',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    },
    "fafb_confident": {
        "files": [
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_confident_v4.csv',
        ],
        "voxel_size": (40, 4, 4),
        "db_name": 'synister_fafb_v4_confident',
        "test_fraction": 0.0,       # fraction of whole dataset to use for testing
        "validation_fraction": 1.0  # fraction of remaining data to use for validation
    },
    "fafb_v3updated": {
        "files": [
            'fafb/consolidated/2021-12-08/2021-11-23_FAFB_connectors_by_hemi_lineage_v3.csv',
        ],
        "holdout_files": [
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_confident_v4.csv',
        ],
        "voxel_size": (40, 4, 4),
        "db_name": 'synister_fafb_v4_v3updated',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    },
    "fafb_v3updatedKC": {
        "files": [
            'fafb/consolidated/2021-12-08/2021-11-23_FAFB_connectors_by_hemi_lineage_v3.csv',
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_v4_KC_only.csv',
        ],
        "holdout_files": [
            'fafb/consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_confident_v4.csv',
        ],
        "voxel_size": (40, 4, 4),
        "db_name": 'synister_fafb_v4_v3updatedKC',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    },
    "hemi": {
        "files": [
            'hemi/consolidated/2021-10-27/hemibrain_connectors_by_hemi_lineage_October2021.json'
        ],
        "voxel_size": (8, 8, 8),
        "db_name": 'synister_hemi_v1',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    },
    "malevnc": {
        "files": [
            'malevnc/consolidated/vnc_filtered_090621/synapses.json'
        ],
        "voxel_size": (8, 8, 8),
        "db_name": 'synister_malevnc_v0',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    }
}

NT_SYNAPSES_THRESHOLD = 1000
NT_SKELETONS_THRESHOLD = 3

parser = argparse.ArgumentParser()
parser.add_argument(
    'dataset',
    type=str,
    help="Name of the dataset to ingest",
    choices=DATASETS.keys())
parser.add_argument(
    '--credentials',
    '-c',
    type=str,
    required=True,
    help="MongoDB credential file")


def read_synapses(synapse_files, voxel_size):

    synapses = []
    for filename in synapse_files:
        with open(filename, 'r') as f:
            synapses += json.load(f)

    # connector_id -> synapse_id

    def get_synapse_id(synapse):

        # use connector_id if available
        if synapse['connector_id'] is not None:
            return synapse['connector_id']

        # fall back to Cantor number of coordinates (int, in voxels)
        return cantor_number(tuple(
            int(synapse[d]) // voxel_size[i]
            for i, d in enumerate(['z', 'y', 'x'])
        ))

    def get_skeleton_id(synapse):
        fields = ["skid", "flywire_id", "body_id"]
        for f in fields:
            val = synapse.get(f)
            if val is not None:
                return val
        return None

    # bring into synapse format as expected by SynisterDb
    synapses = [
        {
            **synapse,
            # synister DB expects int for coordinates
            'x': int(synapse['x']),
            'y': int(synapse['y']),
            'z': int(synapse['z']),
            'synapse_id': get_synapse_id(synapse),
            'skeleton_id': get_skeleton_id(synapse),
            'brain_region': synapse['region']
        }
        for synapse in synapses
    ]

    # filter underrepresented neurotransmitters
    neurotransmitter_counts = defaultdict(lambda: {"synapse": 0, "skeleton": set()})
    for synapse in synapses:
        counts = neurotransmitter_counts[synapse["neurotransmitter"]]
        counts["synapse"] += 1
        counts["skeleton"].add(synapse["skeleton_id"])

    for nt in neurotransmitter_counts.keys():
        neurotransmitter_counts[nt]["skeleton"] = len(neurotransmitter_counts[nt]["skeleton"])
    print("Neurotransmitter counts")
    print(neurotransmitter_counts)

    accepted_neurotransmitters = []
    for nt, c in neurotransmitter_counts.items():
        if c["synapse"] >= NT_SYNAPSES_THRESHOLD and c["skeleton"] >= NT_SKELETONS_THRESHOLD:
            accepted_neurotransmitters.append(nt)
        else:
            print(f"Excluding {nt}")

    filtered_synapses = []
    for synapse in synapses:
        if synapse['neurotransmitter'] not in accepted_neurotransmitters:
            print(f"Skipping {synapse} with filtered neurotransmitter")
            continue
        filtered_synapses.append(synapse)
    synapses = filtered_synapses

    return synapses


def ingest_synapses(synapses, db):

    # check for duplicate IDs
    synapse_ids = np.array([
        synapse['synapse_id']
        for synapse in synapses])
    ids, unique_idxs, inverse_idxs, counts = np.unique(
        synapse_ids, return_index=True, return_inverse=True, return_counts=True)
    multiple = (counts > 1)
    duplicate_ids = ids[multiple]
    duplicate_ids_idxs = unique_idxs[multiple]
    duplicate_counts = counts[multiple]
    nonduplicated_idxs = unique_idxs[~multiple]

    if len(duplicate_ids) > 0:
        print(f"Found {len(duplicate_ids)} duplicate synapse IDs")

        # Find identical duplicates.
        identical_duplicates = np.empty(len(ids), dtype=bool)
        identical_duplicates[multiple] = True
        for uniq_idx, s in zip(inverse_idxs, synapses):
            if not identical_duplicates[uniq_idx]:
                continue
            reference = synapses[unique_idxs[uniq_idx]]
            identical_duplicates[uniq_idx] = reference == s
        identical_duplicates = identical_duplicates[multiple]
        identical_reference_ids = duplicate_ids_idxs[identical_duplicates]

        def print_duplicates(ids, counts):
            for duplicate, count in zip(ids, counts):
                print(f"{duplicate} repeats {count} times:")
                indices = (synapse_ids == duplicate).nonzero()[0]
                for index in indices:
                    print(synapses[index])

        print(f"Of these, {sum(identical_duplicates)} IDs are identical and will be deduplicated:")
        print("(showing at most 100)")
        print_duplicates(
            duplicate_ids[identical_duplicates][:100],
            duplicate_counts[identical_duplicates][:100])

        print(f"The remaining {sum(~identical_duplicates)} IDs are not identical and will be removed:")
        print("(showing at most 100)")
        print_duplicates(
            duplicate_ids[~identical_duplicates][:100],
            duplicate_counts[~identical_duplicates][:100])

        # Deduplicate identical duplicates and remove non-identical duplicates.
        num_deduplicated = sum(duplicate_counts[identical_duplicates]) - len(identical_reference_ids)
        print(f"Removing {num_deduplicated} identical duplicated synapses.")
        print(len(synapses))
        retain_mask = np.zeros(len(synapse_ids), dtype=bool)
        retain_mask[nonduplicated_idxs] = True
        retain_mask[identical_reference_ids] = True
        synapses = np.asarray(synapses)[retain_mask].tolist()
        print(len(synapses))


    # hemi_lineage_id, hemi_lineage_name

    hemi_lineages = {}
    hemi_lineage_id = 0
    for synapse in synapses:
        hemi_lineage_name = synapse['hemilineage']
        if hemi_lineage_name not in hemi_lineages:
            hemi_lineages[hemi_lineage_name] = {
                **db.hemi_lineage,
                'hemi_lineage_name': hemi_lineage_name,
                'hemi_lineage_id': hemi_lineage_id
            }
            hemi_lineage_id += 1
    synister_hemi_lineages = list(hemi_lineages.values())

    # skeleton_id, hemi_lineage_id, nt_known, type=None, match=None, quality=None

    skeletons = {}
    for synapse in synapses:
        skeleton_id = synapse["skeleton_id"]
        if skeleton_id not in skeletons:
            skeletons[skeleton_id] = {
                **db.skeleton,
                'skeleton_id': skeleton_id,
                'hemi_lineage_id': hemi_lineages[synapse['hemilineage']]['hemi_lineage_id'],
                'nt_known': [synapse['neurotransmitter']]
            }
    synister_skeletons = list(skeletons.values())

    # write to DB

    db.write(
        synapses=synapses,
        skeletons=synister_skeletons,
        hemi_lineages=synister_hemi_lineages)


def create_synapse_split(
        synapses,
        split_attribute,
        split_name,
        test_fraction,
        validation_fraction):

    print()
    print()
    print(f"Creating split by {split_attribute}, test fraction = "
          f"{test_fraction}, validation fraction = {validation_fraction}...")

    skipped_attribute = 0
    skipped_nt = 0
    nt_by_synapse_id = {}
    superset_by_synapse_id = {}
    for synapse in synapses:
        if synapse[split_attribute] is None:
            skipped_attribute += 1
            continue
        if synapse['neurotransmitter'] is None:
            skipped_nt += 1
            continue
        synapse_id = synapse['connector_id']
        attribute = synapse[split_attribute]
        neurotransmitter = synapse['neurotransmitter']
        nt_by_synapse_id[synapse_id] = (neurotransmitter,)  # synister wants a tuple
        superset_by_synapse_id[synapse_id] = attribute

    supersets = list(set(superset_by_synapse_id.values()))
    neurotransmitters = list(set(nt_by_synapse_id.values()))
    synapse_ids = list(nt_by_synapse_id.keys())

    print()
    print(f"Found {len(supersets)} different values for {split_attribute}")
    if len(supersets) <= 100:
        print(supersets)
    else:
        print(f"{supersets[:100]} (and {len(supersets) - 100} more...)")
    print(f"Skipped {skipped_attribute}/{len(synapses)} synapses without a "
          f"'{split_attribute}' attribute")
    print(f"Skipped {skipped_nt}/{len(synapses)} synapses without a "
          "'neurotransmitter' attribute")
    print(f"Found neurotransmitters {list([n[0] for n in neurotransmitters])}")

    if not synapse_ids:
        print("No synapses left with the required attributes, skipping split")
        return

    train_validation_synapse_ids, test_synapse_ids, neurotransmitters, synapse_split_nts = find_optimal_split_or_fallback(
        synapse_ids=synapse_ids,
        superset_by_synapse_id=superset_by_synapse_id,
        nt_by_synapse_id=nt_by_synapse_id,
        supersets=supersets,
        neurotransmitters=neurotransmitters,
        synapse_split_nts=[],
        set_b_fraction=test_fraction,
        set_a_name="(train ∪ validation)",
        set_b_name="test",
        split_attribute=split_attribute,
    )

    train_synapse_ids, validation_synapse_ids, neurotransmitters, synapse_split_nts = find_optimal_split_or_fallback(
        synapse_ids=train_validation_synapse_ids,
        superset_by_synapse_id=superset_by_synapse_id,
        nt_by_synapse_id=nt_by_synapse_id,
        supersets=supersets,
        neurotransmitters=neurotransmitters,
        synapse_split_nts=synapse_split_nts,
        set_b_fraction=validation_fraction,
        set_a_name="train",
        set_b_name="validation",
        split_attribute=split_attribute,
    )

    # store split in DB

    db.make_split(
       split_name,
       train_synapse_ids,
       test_synapse_ids,
       validation_synapse_ids)

    return (
        train_synapse_ids,
        test_synapse_ids,
        validation_synapse_ids)


def find_optimal_split_or_fallback(
        synapse_ids,
        superset_by_synapse_id,
        nt_by_synapse_id,
        supersets,
        neurotransmitters,
        synapse_split_nts,
        set_b_fraction,
        set_a_name,
        set_b_name,
        split_attribute,
    ):
        resplit = True

        while resplit:
            resplit = False
            try:

                print()
                print(f"Creating {set_a_name} / {set_b_name} split, "
                    f"objective is {1.0 - set_b_fraction}/{set_b_fraction}:")
                filtered_synapse_ids = [sid for sid in synapse_ids if nt_by_synapse_id[sid] in neurotransmitters]

                if set_b_fraction <= 0.0:

                    print(f"({set_b_name} fraction is 0, no need to split)")
                    a_set_synapse_ids = filtered_synapse_ids
                    b_set_synapse_ids = []

                elif set_b_fraction >= 1.0:

                    print(f"({set_b_name} fraction is 1, no need to split)")
                    a_set_synapse_ids = []
                    b_set_synapse_ids = filtered_synapse_ids

                else:

                    a_set, b_set = find_optimal_split(
                        synapse_ids=filtered_synapse_ids,
                        superset_by_synapse_id=superset_by_synapse_id,
                        nt_by_synapse_id=nt_by_synapse_id,
                        neurotransmitters=neurotransmitters,
                        supersets=supersets,
                        train_fraction=1.0 - set_b_fraction)

                    a_set_synapse_ids = list(itertools.chain(*a_set.values()))
                    b_set_synapse_ids = list(itertools.chain(*b_set.values()))

            except ImpossibleSplit as e:

                print()
                print(
                    f"\tWARNING: failed to create optimal split for {e.nt} on "
                    f"attribute {split_attribute}!")
                print(
                    f"\toptimal fraction {e.optimal_fraction} deviates too far from "
                    f"target {e.target_fraction}")
                print()
                print(f"\tFalling back to {set_a_name} / {set_b_name} split on {e.nt} synapses")
                print()

                resplit = True
                neurotransmitters.remove(e.nt)
                synapse_split_nts.append(e.nt)


        for nt in synapse_split_nts:

            synapse_ids = [
                synapse_id
                for synapse_id, _nt in nt_by_synapse_id.items()
                if _nt == nt
            ]

            random.seed(19120623)
            random.shuffle(synapse_ids)
            split_index = int((1.0 - set_b_fraction) * len(synapse_ids))
            a_set_synapse_ids += synapse_ids[:split_index]
            b_set_synapse_ids += synapse_ids[split_index:]

            print(
                f"Split {nt} randomly per synapse into "
                f"{split_index}/{len(synapse_ids)} = "
                f"{100.0 * split_index/len(synapse_ids)}%")

        return a_set_synapse_ids, b_set_synapse_ids, neurotransmitters, synapse_split_nts

if __name__ == '__main__':

    args = parser.parse_args()
    dataset = DATASETS[args.dataset]

    db = SynisterDb(args.credentials, dataset["db_name"])
    db.create(overwrite=True)

    synapses = read_synapses(dataset['files'], dataset['voxel_size'])

    ingest_synapses(synapses, db)

    has_holdout = "holdout_files" in dataset
    if has_holdout:
        holdout_synapses = read_synapses(dataset['holdout_files'], dataset['voxel_size'])
        holdout_synapse_ids = set(s["synapse_id"] for s in holdout_synapses)
        original_len = len(synapses)
        synapses = [s for s in synapses if not s["synapse_id"] in holdout_synapse_ids]
        print(f"Excluded {original_len - len(synapses)}/{original_len} holdout synapses.")

    db.init_splits()

    create_synapse_split(
        synapses,
        'skeleton_id',
        'skeleton',
        test_fraction=dataset['test_fraction'],
        validation_fraction=dataset['validation_fraction'])

    snt_splits = create_synapse_split(
        synapses,
        'skeleton_id',
        'skeleton_no_test',
        test_fraction=0.0,
        validation_fraction=dataset['validation_fraction'])

    if has_holdout:
        snt_splits[0].extend(holdout_synapse_ids)
        db.make_split(
            "skeleton_no_test_including_holdout",
            *snt_splits,
        )

    create_synapse_split(
        synapses,
        'brain_region',
        'brain_region',
        test_fraction=dataset['test_fraction'],
        validation_fraction=dataset['validation_fraction'])
