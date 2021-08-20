from configparser import ConfigParser
from funlib.math import cantor_number
from synister import SynisterDb, find_optimal_split
import argparse
import json
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    '--credentials',
    '-c',
    type=str,
    required=True,
    help="MongoDB credential file")

if __name__ == '__main__':

    files = [
        'fafb/consolidated/2021-06-11/FAFB_connectors_by_hemi_lineage_June2021.json',
        'fafb/consolidated/2021-06-11/FAFB_connectors_with_known_neurotransmitter_by_compartment_2020.json',
        'fafb/consolidated/2021-06-11/FAFB_verified_predicted_synapses_by_transmitter_June2021.json'
    ]
    voxel_size = (40, 4, 4)
    db_name = 'synister_fafb_v4'

    args = parser.parse_args()
    db = SynisterDb(args.credentials, db_name)
    db.create(overwrite=True)


    synapses = []
    for filename in files:
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

    synister_synapses = [
        {
            # synister DB expects int for coordinates
            'x': int(synapse['x']),
            'y': int(synapse['y']),
            'z': int(synapse['z']),
            'synapse_id': get_synapse_id(synapse),
            'skeleton_id': synapse.get('skid', synapse.get('flywire_id', synapse.get('body_id', None))),
            'compartment': synapse['compartment'],
            'brain_region': [synapse['region']]
        }
        for synapse in synapses
    ]

    # check for duplicate IDs
    synapse_ids = np.array([
        synapse['synapse_id']
        for synapse in synister_synapses])
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
                'skeleton_id': skeleton_id,
                'hemi_lineage_id': hemi_lineages[synapse['hemilineage']]['hemi_lineage_id'],
                'nt_known': [synapse['neurotransmitter']]
            }
    synister_skeletons = list(skeletons.values())

    # write to DB

    db.write(
        synapses=synister_synapses,
        skeletons=synister_skeletons,
        hemi_lineages=synister_hemi_lineages)

    # create split by brain region


    neurotransmitters = [
        ('gaba',),
        ('acetylcholine',),
        ('glutamate',),
        ('dopamine',),
        ('octopamine',),
        ('serotonin',),
    ]

    region_by_synapse_id = {}
    nt_by_synapse_id = {}
    synapse_ids = []

    regions = list(set([
        synapse['region']
        for synapse in synapses
    ]))

    print(f"Creating split by region, for regions: {regions}")

    skipped_region = 0
    skipped_nt = 0
    for synapse in synapses:
        if synapse['region'] is None:
            skipped_region += 1
            continue
        if (synapse['neurotransmitter'],) not in neurotransmitters:
            skipped_nt += 1
            continue
        region_by_synapse_id[synapse['connector_id']] = synapse['region']
        nt_by_synapse_id[synapse['connector_id']] = (synapse['neurotransmitter'],)
        synapse_ids.append(synapse['connector_id'])

    print(f"Skipped {skipped_region}/{len(synapses)} synapses without a region attribute")
    print(f"Skipped {skipped_nt}/{len(synapses)} synapses without a NT in {neurotransmitters}")

    # split into train and test

    train_set, test_set = find_optimal_split(
        synapse_ids=synapse_ids,
        superset_by_synapse_id=region_by_synapse_id,
        nt_by_synapse_id=nt_by_synapse_id,
        neurotransmitters=neurotransmitters,
        supersets=regions,
        train_fraction=0.8)

    # split train further into train and validate

    train_synapse_ids = []
    for ids in train_set.values():
        train_synapse_ids += ids

    train_set, validation_set = find_optimal_split(
        synapse_ids=train_synapse_ids,
        superset_by_synapse_id=region_by_synapse_id,
        nt_by_synapse_id=nt_by_synapse_id,
        neurotransmitters=neurotransmitters,
        supersets=regions,
        train_fraction=0.8)

    # convert train, validate, and test sets into lists of synapse IDs

    train_synapse_ids = []
    for ids in train_set.values():
        train_synapse_ids += ids
    validation_synapse_ids = []
    for ids in validation_set.values():
        validation_synapse_ids += ids
    test_synapse_ids = []
    for ids in test_set.values():
        test_synapse_ids += ids

    # store split in DB

    db.make_split(
       'brain_region',
       train_synapse_ids,
       test_synapse_ids,
       validation_synapse_ids)
