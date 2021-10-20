import numpy as np
import copy
from synister.synister_db import SynisterDb

def ingest(db_credentials, db_name, synapses, skeletons, overwrite=False):
    db = SynisterDb(db_credentials, db_name)
    db.create(overwrite=overwrite)
    db.write(synapses, skeletons)

def read_csv(file_path, skeleton_id_column, synapse_id_column, x_column, y_column, z_column, nt=None, voxel_size=np.array([8,8,8])):
    """
    Reads specified columns from a csv and returns list of synapses and skeleton dicts 
    ready to be ingested into the database.
    """
    name_map = {synapse_id_column: "synapse_id", skeleton_id_column: "skeleton_id", x_column: "x", y_column: "y", z_column: "z"}
    data = np.genfromtxt(file_path, delimiter=",", dtype=None, encoding="utf-8", case_sensitive="lower")
    header = list(data[0])

    name_to_column = {"synapse_id": header.index(synapse_id_column), 
                      "skeleton_id": header.index(skeleton_id_column),
                      "x": header.index(x_column),
                      "y": header.index(y_column),
                      "z": header.index(z_column)}

    synapse = {"x": None, "y": None, "z": None, "synapse_id": None, "skeleton_id": None}
    skeleton = {"skeleton_id": None, "nt_known": None, "hemi_lineage_id": None}


    synapses = []
    skeletons = []
    skids = set([])

    for row in data[1:]:
        synapse_tmp = copy.copy(synapse)
        for key in synapse:
            synapse_tmp[key] = int(row[name_to_column[key]])

        synapse_tmp = make_physical(voxel_size,synapse_tmp)
        synapses.append(synapse_tmp)
        
        skid = int(row[name_to_column["skeleton_id"]])
        if not skid in skids:
            skids.add(skid)
            skeleton_tmp = copy.copy(skeleton)
            skeleton_tmp["skeleton_id"] = skid
            skeleton_tmp["nt_known"] = [nt]
            skeletons.append(skeleton_tmp)

    return synapses, skeletons

def make_physical(voxel_size, synapse):
    coords = ["x", "y", "z"]

    for k in range(3):
        synapse[coords[k]] = int(synapse[coords[k]] * voxel_size[2-k])

    return synapse

def write_vnc():
    for nt in ["gaba", "acetylcholine", "glutamate"]:
        synapses, skeletons = read_csv(f"/nrs/funke/ecksteinn/nils_data/synister_data/database_source_data/vnc_filtered_090621/{nt}.csv",
                                       "body", "syn_id", "x", "y", "z", nt=nt)
        print(len(synapses), len(skeletons))

        ingest(db_credentials="/nrs/funke/ecksteinn/nils_data/synister_data/credentials/db_credentials.ini", 
               db_name="synister_vnc_v0", 
               synapses=synapses, 
               skeletons=skeletons, 
               overwrite=False)

    
if __name__ == "__main__":
    write_vnc()
