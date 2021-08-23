from csv import DictReader
import json

hemi_lineage_in_file = 'original/2021-08-21/FAFB_connectors_by_hemi_lineage_August2021.csv'
hemi_lineage_out_file = 'consolidated/2021-08-21/FAFB_connectors_by_hemi_lineage_August2021.json'

verified_in_file = 'original/2021-08-21/FAFB_verified_predicted_synapses_by_transmitter_August2021.csv'
verified_out_file = 'consolidated/2021-08-21/FAFB_verified_predicted_synapses_by_transmitter_August2021.json'

def read_csv(
        filename,
        classic_other=False,
        verified_column=False,
        flywire=False):

    print(f"Reading {filename}")

    unexpected_nt_types = set()

    with open(filename, 'r') as f:

        reader = DictReader(f)

        skids = set()
        skip_skids = set()

        synapses = []

        for row in reader:

            if flywire:
                skid = int(row['flywire.id'])
            else:
                skid = int(row['skid'])
            skids.add(skid)

            if skid in skip_skids:
                continue

            if classic_other:

                classic = row['known.classic.transmitter']
                other = row['known.other.transmitter']

                neurotransmitter = (
                    classic if classic != 'unknown' else (
                        other if other != 'unknown' else None
                    )
                )
            else:
                if 'transmitter' in row:
                    neurotransmitter = row['transmitter']
                else:
                    neurotransmitter = row['known.neurotransmitter']
                if neurotransmitter == 'unknown':
                    neurotransmitter = None

            if neurotransmitter is not None:
                neurotransmitter = neurotransmitter.lower()
                if neurotransmitter == 'kc_acetylcholine':
                    neurotransmitter = 'acetylcholine'

            if verified_column:
                verified = (row['neurotransmitter.verified'].lower() == 'true')
                if not verified:
                    neurotransmitter = None

            if neurotransmitter is None:
                skip_skids.add(skid)
                print(f"Skipping skeleton {skid} (no known NT)")
                continue

            if neurotransmitter not in [
                    'gaba',
                    'acetylcholine',
                    'glutamate',
                    'octopamine',
                    'serotonin',
                    'dopamine']:
                unexpected_nt_types.add(neurotransmitter)

            hemi_lineage = row['ItoLee.Hemilineage']
            lineage = row['ItoLee.Lineage']
            if hemi_lineage == '' or hemi_lineage == 'NA' or hemi_lineage == 'unknown':
                hemi_lineage = None
            if lineage == '' or lineage == 'NA' or lineage == 'unknown':
                lineage = None

            compartment = None
            if 'compartment' in row:
                compartment = row['compartment']

            connector_id = row['connector_id']
            if connector_id == 'unknown':
                connector_id = None
            else:
                connector_id = int(connector_id)

            region = None
            if 'inside' in row:
                region = row['inside']
                if region == 'unknown':
                    region = None

            try:
                x = float(row['x'])
                y = float(row['y'])
                z = float(row['z']) - 40.0  # from CATMAID v14 to FAFB v14 N5 space
            except ValueError as e:
                print(
                    f"Error parsing coordinates of synapse {connector_id} "
                    f"(skid = {skid}): {e}")
                print("Skipping this synapse")
                continue

            flywire_id = None
            if flywire:
                flywire_id = skid
                skid = None

            synapses.append({
                'skid': skid,
                'flywire_id': flywire_id,
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

        print(f"Skipped {len(skip_skids)}/{len(skids)} skeletons")

        print(f"Encountered unexpected NT types: {unexpected_nt_types}")
        return synapses


synapses = read_csv(hemi_lineage_in_file, classic_other=True)
with open(hemi_lineage_out_file, 'w') as f:
    json.dump(synapses, f, indent=2)

synapses = read_csv(verified_in_file, flywire=True)
with open(verified_out_file, 'w') as f:
    json.dump(synapses, f, indent=2)
