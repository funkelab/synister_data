from csv import DictReader
import json

in_file = 'original/hemibrain_connectors_by_hemi_lineage_June2021_filtered.csv'
out_file = 'consolidated/hemibrain_connectors_by_hemi_lineage_June2021_filtered.json'

def read_csv(filename):

    print(f"Reading {filename}")

    with open(filename, 'r') as f:

        reader = DictReader(f)

        body_ids = set()
        skip_body_ids = set()

        synapses = []

        for row in reader:

            body_id = int(row['bodyid'])
            body_ids.add(body_id)

            if body_id in skip_body_ids:
                continue


            classic = row['known.classic.transmitter']
            other = row['known.other.transmitter']

            neurotransmitter = (
                classic if classic != 'unknown' else (
                    other if other != 'unknown' else None
                )
            )

            if neurotransmitter is not None:
                neurotransmitter = neurotransmitter.lower()

            if neurotransmitter is None:
                skip_body_ids.add(body_id)
                print(f"Skipping skeleton {body_id} (no known NT)")
                continue

            hemi_lineage = row['ItoLee.Hemilineage']
            lineage = None
            if hemi_lineage == '' or hemi_lineage == 'NA' or hemi_lineage == 'unknown':
                hemi_lineage = None

            compartment = None
            if 'Label' in row:
                compartment = row['Label']

            connector_id = row['connector_id']
            if connector_id == 'unknown':
                connector_id = None
            else:
                connector_id = int(connector_id)

            region = row['inside']

            try:
                x = float(row['x'])
                y = float(row['y'])
                z = float(row['z'])
            except ValueError as e:
                print(
                    f"Error parsing coordinates of synapse {connector_id} "
                    f"(body_id = {body_id}): {e}")
                print("Skipping this synapse")
                continue

            synapses.append({
                'body_id': body_id,
                'connector_id': connector_id,
                'x': x,
                'y': y,
                'z': z,
                'hemilineage': hemi_lineage,
                'lineage': lineage,
                'compartment': compartment,
                'region': region,
                'neurotransmitter': neurotransmitter
            })

        print(f"Skipped {len(skip_body_ids)}/{len(body_ids)} skeletons")
        return synapses


synapses = read_csv(in_file)
with open(out_file, 'w') as f:
    json.dump(synapses, f, indent=2)
