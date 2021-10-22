from configparser import ConfigParser
from funlib.math import cantor_number
from synister import SynisterDb, find_optimal_split
import argparse
import json
import numpy as np

DATASETS = {
    "fafb": {
        "files": [
            'fafb/consolidated/2021-08-21/FAFB_connectors_by_hemi_lineage_August2021.json',
            'fafb/consolidated/2021-08-21/FAFB_verified_predicted_synapses_by_transmitter_August2021.json'
        ],
        "voxel_size": (40, 4, 4),
        "db_name": 'synister_fafb_v4',
        "test_fraction": 0.2,       # fraction of whole dataset to use for testing
        "validation_fraction": 0.2  # fraction of remaining data to use for validation
    },
    "hemi": {
        "files": [
            'hemi/consolidated/2021-08-21/hemibrain_connectors_by_hemi_lineage_August2021.json'
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

    # bring into synapse format as expected by SynisterDb
    synapses = [
        {
            **synapse,
            # synister DB expects int for coordinates
            'x': int(synapse['x']),
            'y': int(synapse['y']),
            'z': int(synapse['z']),
            'synapse_id': get_synapse_id(synapse),
            'skeleton_id': synapse.get(
                'skid',
                synapse.get(
                    'flywire_id',
                    synapse.get('body_id', None))),
            'brain_region': synapse['region']
        }
        for synapse in synapses
    ]

    return synapses


def ingest_synapses(synapses, db):

    # check for duplicate IDs
    synapse_ids = np.array([
        synapse['synapse_id']
        for synapse in synapses])
    ids, counts = np.unique(synapse_ids, return_counts=True)
    multiple = (counts > 1)
    duplicate_ids = ids[multiple]
    duplicate_counts = counts[multiple]
    if len(duplicate_ids) > 0:
        print(f"Found {len(duplicate_ids)} duplicate synapse IDs")
        print("(showing at most 100)")
        for duplicate, count in zip(
                duplicate_ids[:100],
                duplicate_counts[:100]):
            print(f"{duplicate} repeats {count} times:")
            same_mask = (synapse_ids == duplicate)
            indices = np.arange(0, len(synapses))[same_mask]
            for index in indices:
                print(synapses[index])


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
        skeleton_id = synapse.get('skid', synapse.get('flywire_id', synapse.get('body_id', None)))
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

    # split off test fraction
    print()
    print("Creating (train ∪ validation) / test split, "
          f"objective is {1.0 - test_fraction}/{test_fraction}:")

    if test_fraction > 0.0:

        train_validation_set, test_set = find_optimal_split(
            synapse_ids=synapse_ids,
            superset_by_synapse_id=superset_by_synapse_id,
            nt_by_synapse_id=nt_by_synapse_id,
            neurotransmitters=neurotransmitters,
            supersets=supersets,
            train_fraction=1.0 - test_fraction)

        train_validation_synapse_ids = []
        for ids in train_validation_set.values():
            train_validation_synapse_ids += ids
        test_synapse_ids = []
        for ids in test_set.values():
            test_synapse_ids += ids

    else:

        print("(test fraction is 0, no need to split)")
        train_validation_synapse_ids = synapse_ids
        test_synapse_ids = []

    # split train_validation_set into train and validation
    print()
    print("Creating train / validation split, "
          f"objective is {1.0 - validation_fraction}/{validation_fraction}:")
    train_set, validation_set = find_optimal_split(
        synapse_ids=train_validation_synapse_ids,
        superset_by_synapse_id=superset_by_synapse_id,
        nt_by_synapse_id=nt_by_synapse_id,
        neurotransmitters=neurotransmitters,
        supersets=supersets,
        train_fraction=1.0 - validation_fraction)

    # convert train, validation, and test sets into lists of synapse IDs

    train_synapse_ids = []
    for ids in train_set.values():
        train_synapse_ids += ids
    validation_synapse_ids = []
    for ids in validation_set.values():
        validation_synapse_ids += ids

    # store split in DB

    db.make_split(
       split_name,
       train_synapse_ids,
       test_synapse_ids,
       validation_synapse_ids)

if __name__ == '__main__':

    args = parser.parse_args()
    dataset = DATASETS[args.dataset]

    db = SynisterDb(args.credentials, dataset["db_name"])
    db.create(overwrite=True)

    synapses = read_synapses(dataset['files'], dataset['voxel_size'])
    ingest_synapses(synapses, db)

    db.init_splits()

    create_synapse_split(
        synapses,
        'skeleton_id',
        'skeleton',
        test_fraction=dataset['test_fraction'],
        validation_fraction=dataset['validation_fraction'])

    create_synapse_split(
        synapses,
        'skeleton_id',
        'skeleton_no_test',
        test_fraction=0.0,
        validation_fraction=dataset['validation_fraction'])

    create_synapse_split(
        synapses,
        'brain_region',
        'brain_region',
        test_fraction=dataset['test_fraction'],
        validation_fraction=dataset['validation_fraction'])
