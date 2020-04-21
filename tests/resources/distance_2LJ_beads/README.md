# Info
Short calculation with LAMMPS for two Lennard-Jones beads pulled apart

* System : 
  * 2 atoms in a box
  * `box block 0 30 0 30 0 30`
* Simulation :
  * Interactions:
    * `units lj`
    * `timestep 0.005`
    * `mass 1 1.0`
    * `pair_style lj/cut 2.5`
    * `pair_coeff 1 1 1.0 1.0 2.5`
  * no integrator, no velocities
  * Shift `atom2` from `min_dist` to `boxsize-min_dist`
    * `step_size` : `0.1`
    * `min_dist` :  `0.4`
    * `x : 0.4 -> 29.6`
    * calculate the potential energy
      * store the data in `dump.lammpstrj.gz`

## Files
* `input.distance_2LJ_beads.lmp` : LAMMPS input file
* `log.lammps` : LAMMPS logfile 
  * properties / step (`Step Temp PotEng TotEng Press Volume 
   E_bond E_angle E_dihed E_impro E_vdwl E_coul E_long E_tail v_e2body v_ecoul_tot`)
    * `e2body    equal ebond+evdwl`
    * `ecoul_tot equal ecoul+elong`
* `dump.lammpstrj` : LAMMPS dump file
  * properties / atom(`id type x y z fx fy fz c_peatom`)

## HowTo
run lammps with:
```bash
lmp -in input.distance_2LJ_beads.lmp
```

* LAMMPS : `stable_3Mar2020` 