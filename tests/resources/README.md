# Info
Resources for the test
* Simulation performed with LAMMPS (https://lammps.sandia.gov/)
    * `distance_2LJ_beads` : Short calculation with LAMMPS for two Lennard-Jones beads pulled apart
    * `LJ_melt` : Short calculation with LAMMPS for a Lennard Jones melt. 


## Folders

* `distance_2LJ_beads` <br>
  Short calculation with LAMMPS for two Lennard-Jones beads pulled apart
  * System : 2 atoms in a box
  * Simulation : 
    * Shift `atom2` from `min_dist` to `boxsize-min_dist`
    * Lennard-Jones interactions & reduced units 

* `LJ_melt` <br>
  Short calculation with LAMMPS for a Lennard Jones melt. 
  * System : 4000 atoms in a box
  * Simulation : 
   * Lennard-Jones interactions & reduced units 