# Set-Up mmcDNSFoam Case

In the constant folder all files for the DNS and LES settings are included. 
The LES settings are found in the filterMesh directory, which also contains
the cloudProperties file for the mmc stochastic particles. 

As mmc particles and Eulerian field should use the same chemistry and thermophysical
models, a sym-link to the thermophysicalProperties and chemistryProperties
file makes sense. 

