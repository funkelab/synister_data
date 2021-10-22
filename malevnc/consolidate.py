from csv import DictReader
import json

ach_in_file = 'original/vnc_filtered_090621/acetylcholine.csv'
gaba_in_file = 'original/vnc_filtered_090621/gaba.csv'
glut_in_file = 'original/vnc_filtered_090621/glutamate.csv'

out_file = 'consolidated/vnc_filtered_090621/synapses.json'

def read_csv(
        filename,
        neurotransmitter,
        verified_column=False,
        flywire=False):

    print(f"Reading {filename}")

    unexpected_nt_types = set()

    with open(filename, 'r') as f:

        reader = DictReader(f)
        synapses = []

        for row in reader:

            skeleton_id = int(row['body'])
            synapse_id = int(row['syn_id'])

            if neurotransmitter not in [
                    'gaba',
                    'acetylcholine',
                    'glutamate',
                    'octopamine',
                    'serotonin',
                    'dopamine']:
                unexpected_nt_types.add(neurotransmitter)

            # we don't have this information
            hemi_lineage = None
            lineage = None
            compartment = None
            region = None

            try:
                # convert from voxels to nm
                x = float(row['x']) * 8
                y = float(row['y']) * 8
                z = float(row['z']) * 8
            except ValueError as e:
                print(
                    f"Error parsing coordinates of synapse {synapse_id} "
                    f"(skeleton_id = {skeleton_id}): {e}")
                print("Skipping this synapse")
                continue

            synapses.append({
                'body_id': skeleton_id,
                'connector_id': synapse_id,
                'x': x,
                'y': y,
                'z': z,
                'hemilineage': hemi_lineage,
                'lineage': lineage,
                'compartment': compartment,
                'region': region,
                'neurotransmitter': neurotransmitter
            })

        print(f"Encountered unexpected NT types: {unexpected_nt_types}")
        return synapses


synapses = read_csv(ach_in_file, neurotransmitter='acetylcholine')
synapses += read_csv(gaba_in_file, neurotransmitter='gaba')
synapses += read_csv(glut_in_file, neurotransmitter='glutamate')

with open(out_file, 'w') as f:
    json.dump(synapses, f, indent=2)
