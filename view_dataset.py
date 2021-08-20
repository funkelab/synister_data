from funlib.show.neuroglancer import add_layer
import configargparse
import daisy
import json
import neuroglancer

parser = configargparse.ArgParser()
parser.add(
    '--raw-container',
    help="Path to the raw data container (zarr or N5)",
    required=True)
parser.add(
    '--raw-dataset',
    help="Name of the raw dataset in the container",
    required=True)
parser.add(
    '--synapse-dataset',
    help="JSON file of synapses to show",
    required=True)

if __name__ == '__main__':

    options = parser.parse_args()

    raw = daisy.open_ds(
        options.raw_container,
        options.raw_dataset)
    print(f"Found raw data in roi {raw.roi}, voxel size {raw.voxel_size}")

    with open(options.synapse_dataset, 'r') as f:
        synapses = json.load(f)
    synapse_annotations = [
        neuroglancer.PointAnnotation(
            id=i,
            point=[
                synapses[i]['z'],
                synapses[i]['y'],
                synapses[i]['x']
            ])
        for i in range(100)
    ]
    print(f"Showing {len(synapse_annotations)} synapses in neuroglancer")
    print("Positions (zyx, in nm) of first 10 synapses:")
    for i in range(10):
        print(list(synapses[i][k] for k in ['z', 'y', 'x']))

    neuroglancer.set_server_bind_address('0.0.0.0')
    viewer = neuroglancer.Viewer()
    with viewer.txn() as s:
        add_layer(s, raw, 'raw')
        dimensions = neuroglancer.CoordinateSpace(
            names=['z', 'y', 'x'],
            units='nm',
            scales=(1, 1, 1))
        s.layers.append(
            name="synapses",
            layer=neuroglancer.LocalAnnotationLayer(
                dimensions=dimensions,
                annotations=synapse_annotations))

    print(viewer)
    input("Press ENTER to quit")
