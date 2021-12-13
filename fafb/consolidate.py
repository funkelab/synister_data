from csv import DictReader
import json
import os

out_path = "consolidated/2021-12-08"
files = {
    "confident": {
        "in_file": "original/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_confident_v4.csv",
        "kwargs": {"classic_other": True},
    },
    "v3_updated": {
        "in_file": "original/2021-11-23/2021-11-23_FAFB_connectors_by_hemi_lineage_v3.csv",
        "kwargs": {"classic_other": True},
    },
    "catmaid": {
        "in_file": "original/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_v4.csv",
        "kwargs": {"classic_other": True},
    },
    "catmaid_kcs_only": {
        "in_file": "original/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_v4.csv",
        "out_file": "consolidated/2021-12-08/2021-12-08_FAFB_connectors_by_hemi_lineage_v4_KC_only.csv",
        "kwargs": {"classic_other": True, "kc_only": True},
    },
    "flywire": {
        "in_file": "original/2021-12-08/2021-12-08_FAFB_verified_predicted_synapses_by_transmitter_v4.csv",
        "kwargs": {"flywire": True},
    },
}

def read_csv(
        filename,
        classic_other=False,
        putative_other=False,
        verified_column=False,
        flywire=False,
        kc_only=False):

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

                if neurotransmitter is None and putative_other:

                    classic = row['putative.classic.transmitter']
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
                elif kc_only:
                    continue

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


for file_description, file in files.items():
    print(f"Consolidating {file_description} synapses...")
    synapses = read_csv(file["in_file"], **file["kwargs"])

    if "out_file" in file:
        out_file = file["out_file"]
    else:
        base_file = os.path.basename(file["in_file"])
        out_file = os.path.join(out_path, base_file)

    with open(out_file, 'w') as f:
        json.dump(synapses, f, indent=2)
