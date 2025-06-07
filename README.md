**<font size=5> A tool to export XYZ trajectory file from Material Studio</font>**




# Usage

* Linux

  ```
  chmod a+x MS_dump
  ./MS_dump system.pdb system.trj
  ```

* Windows

  ```
  .\MS_dump.exe system.pdb system.trj
  ```



Output `MS_traj.xyz`, the comment line has box information that can be read by `Ovito` directly. 

> Lattice="14.408798 0.0 0.0 0.000000 14.408798 0.0 0.000000 0.000000 14.408798" Properties=species:S:1:pos:R:3



**Note:**

> The exported xyz will use `C` name for all atoms if not provide pdb file, such as:
>
> ```
> .\MS_dump.exe system.trj
> ```

