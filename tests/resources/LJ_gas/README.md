# Info
Short calculation with LAMMPS for a Lennard Jones gas. 
Originally `lammps/examples/melt/in.melt`.

* System : 
  * 32 atoms in a box
  * box : `(6.839903786706787, 6.839903786706787, 6.839903786706787, 90, 90, 90)`
* Simulation :
  * Interactions:
    * `units lj`
    * `timestep 0.005`
    * `mass 1 1.0`
    * `pair_style lj/cut 2.5`
    * `pair_coeff 1 1 1.0 1.0 2.5`
  * Integrator : Velocity Verlet
  * Temperature : `0.3`
  * store the data in `dump.lammpstrj.gz`

## Files
* `input.gas.lmp` : LAMMPS input file
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
lmp -in input.gas.lmp
```

* LAMMPS : `stable_3Mar2020` 