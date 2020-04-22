# Info
Resources for the test
* Simulation performed with LAMMPS (https://lammps.sandia.gov/)
    * `distance_2LJ_beads` : Short calculation with LAMMPS for two Lennard-Jones beads pulled apart
    * `LJ_melt` : Short calculation with LAMMPS for a Lennard Jones melt. 
    * `LJ_gas` : Short calculation with LAMMPS for a Lennard Jones gas. 

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
   
 * `LJ_gas` <br>
  Short calculation with LAMMPS for a Lennard Jones gas. 
  * System : 32 atoms in a box
  * Simulation : 
   * Lennard-Jones interactions & reduced units 