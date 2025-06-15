**<font size=5> A tool to export common trajectory file from Materials Studio</font>**



# Features

* Support full periodic boundary conditions

  > Please note that this tool has converted the coordinates to match the PDB unit cell and coordinates exported by MS, so there may be differences from the original data output by Perl/Trj2Ascii, especially for the triclinic system.

* Support export move + fix atoms if exist `.pdb`

* Support export `xyz` trajectory file

* Support export `xtc` of gromacs file (includes time and step)

* Support export `.trr` format which contains velocities (`nm/ps`) and forces`(kJ/mol/nm)` if it exists

  > Note that fixed atoms will have zero velocity and force, meaning their coordinates/velocity/force data won't appear in the `.trj` file.

* Support export `.txt` format (`Plain text file`) which contains all kinds of energy items and others:

  ```
  #Time(ps) Temperature(K) Potential(kJ/mol) Kinetic(KJ/mol) TotalEnergy(KJ/mol) Pressure(bar) Volume(A^3)
  ```
* Support .whl for Python API
  ```
  pip install xxx.whl
  ```
  Usage reference: [api_test](https://github.com/liuyujie714/MS_Trajdump/blob/master/PyMSDump/api_test.py)
  

# Usage

First locate molecular dynamics trajectory file created by `Materials Studio`, hidden file `.trj`  and `.xtd` are located in same folder.



Then download program from here: [Download](https://github.com/liuyujie714/MS_Trajdump/releases)




* Linux

  ```
  chmod a+x MS_dump
  ./MS_dump -s system.pdb -f system.trj
  ```

* Windows

  ```
  .\MS_dump.exe -s system.pdb -f system.trj
  ```



Default output `MS_traj.xyz`, the comment line has box information that can be read by [Ovito](https://www.ovito.org/)  software directly. 

> Lattice="14.408798 0.0 0.0 0.000000 14.408798 0.0 0.000000 0.000000 14.408798" Properties=species:S:1:pos:R:3





**Note:**

> The exported xyz will use `C` name for all atoms if not provide pdb file, such as:
>
> ```
> .\MS_dump.exe -f system.trj
> ```





`-o` option can control output format



If you want to export `.xtc/.trr` of gromacs, use command:

```
.\MS_dump.exe -s system.pdb -f system.trj -o system.xtc
```

```
.\MS_dump.exe -s system.pdb -f system.trj -o system.trr
```



Also export energy items:

```
.\MS_dump.exe -f system.trj -o energy.txt
```

