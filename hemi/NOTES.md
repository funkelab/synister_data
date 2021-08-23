Original data file from https://www.dropbox.com/s/xghzodj01stxmb7/hemibrain_connectors_by_hemi_lineage_June2021.csv?dl=0

Alex: Here are the new hemibrain synapse GT. This would replace previous files
I shared with Nils for this purpose. Transmitter assignments are based on
manually matching neurons to our CATMAID GT. As with the last file, there are a
lot of unknown fields by design that you’ll end up ignoring. Each row is a
synapse. The bodyid column = hemibrain neuron ID, the skid column = CATMAID
neuron ID. Synapses in hemibrain raw voxel space.

2021/08/21
==========

Update from Alex:

1. Missing hemilineages in Hemibrain data ---> In the FAFB dataset, 12% of
   synapses have an unknown hemilineage. In the Hemibrain one, 23% are missing
   for neurons that were founds by matching to CATMAID data. However, since I
   was updating the GT today, I did a little more research and added even more
   hemibrain synapses in some of our underepresented bins + a few more unusual
   transmitters (e.g. NPF, Drosulfakinin, glycin) in case anything can be done
   with them. This pushes the missing hemilineage count to 47%.  There are now
   5144573 hemibrain synapses. As I posted above, not sure what you want to do
   with the range of confidence scores, but all these synapses or on
   idetnifiable axons/dendrites so I have already cut out a lot of chaff, i.e.
   stuff on somas, primary neurites, outside the neuropil volume, etc.

New data in 2021-08-21 directory

• one file: `hemibrain_connectors_by_hemi_lineage_August2021.csv`
  ⇒ no need to filter any longer
