# suffosion
Numerical simulation of suffosion in porous media

This is numerical study of a coupled erosion and deposition process. 
The mathematical model is based on the theory of flow in deformable porous media with mass-variable skeleton. 
3-D numerical simulation of water injection into the layer is performed. 
It is shown that precipitation of particles which were eroded from solid phase by fluid flow near a well (suffosion) 
causes clogging of porous space at some distance from the injecting well where flow velocity is lower.
Parallel AMGCL solver library was used for calculations on a GPU. 
New block matrix preconditioner was developed in order fast solution of the coupled displacements-pressure linear system.
