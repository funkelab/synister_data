2021/08/13
==========

Source Files
------------

## `FAFB_connectors_by_hemi_lineage_June2021`

• 2021 version supersedes the 2020 version
  ⇒ Alex: can ignore 2020 version
    ! however, 2021 contains less synapses than 2020
• known.classic and known.other
• 204946 synpases without known NT (both classic and other 'unknown')
• coordinates in nm

## `FAFB_connectors_with_known_neurotransmitter_by_compartment_2020`

• known.neurotransmitter and neurotransmitter.verified (a bool)
  ⇒ ignore all where verified = FALSE (those have ACh and GABA as types)
• coordinates in nm

## `FAFB_verified_predicted_synapses_by_transmitter_June2021`

• "transmitter" column
• comes from FlyWire neurons
• coordinates are in CATMAID v14 space, in voxels
  ⇒ z is in sections, presumably needs to be z corrected for N5
  ⇒ some x/y are "unknown"

For All Coordinates
-------------------

• need to subtract 40nm from z to convert from CATMAID to FAFB N5
  ⇒ unclear what to do for FlyWire coordinates?

2021/08/20
==========

Manually check locations of consolidated synapse JSONs

✔ `2021-06-11/FAFB_connectors_by_hemi_lineage_June2021.json`
  ⇒ checked first ~100 synapses
  ✔ annotations z, y, x in nm
  ✔ annotations close to t-bars
  ✔ annotations clearly on pre-synaptic site

✔ `FAFB_connectors_with_known_neurotransmitter_by_compartment_2020.json`
  ⇒ checked first ~100 synapses
  ✔ annotations z, y, x in nm
  ✔ annotations close to t-bars
  ✔ annotations clearly on pre-synaptic site

✔ `FAFB_verified_predicted_synapses_by_transmitter_June2021.json`
  ⇒ `flywire_id` is not a synapse ID, but a neuron ID
    ✔ fixed ingest script
  ⇒ -40nm offset okay here?
    ✔ manually checked a few synapses, seems okay with the -40nm offset
      (supposedly the v14 coordinates have been computed for the CATMAD
      reference)
  ✔ annotations z, y, x in nm
  ✔ annotations close to t-bars
  ✔ annotations clearly on pre-synaptic site

2021/08/21
==========

Update from Alex:

1. Problem with missing synapses versus 2020 ---> The issue was that some of my
   CATMAID reads failed, and I forgot to try to re-fetch neurons that had
   failed. Corrected. Our new corpus is 465757 synapses across 5043 neurons.
   The 2020 file had 339506 synapses across 3019 neurons. This was a critical
   issue, thanks for catching @Jan Funke.

2. Missing brain region for some synapses ---> I recalculated and improved the
   number that gets assigned to a neuropil in the FAFB data. Turns out a
   proportion of synapses cannot be assigned easily (see image, though it’s
   tough in 3D) because they fall just outside meshes or between them. We could
   assign them to the correct mesh by nearest neighbour but it’s only
   7238/465757 synapses, and so I think better to just drop them from brain
   region comparison analyses? (The region names between FAFB and the hemibrain
   data are not exactly the same, because things were named slightly
   differently in the hemibrain native meshes, but they both obey the Ito et
   al. 2014 naming conventions so are mappable.)

New data in 2021-08-21 directory

Only two files now:

  • `FAFB_connectors_by_hemi_lineage_August2021.csv`
    ⇒ contains synapses from `FAFB_connectors_with_known_neurotransmitter_by_compartment_2020.csv`
      (the latter file is not needed any more)
  • `FAFB_verified_predicted_synapses_by_transmitter_August2021.csv`
