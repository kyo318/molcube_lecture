
# 1. 환경설정

[polyply 공식 링크](https://github.com/marrink-lab/polyply_1.0)

```bash
conda create -n test python=3.9 -y 

conda activate test

pip install polyply

(pip install numba)
```

- 사용 가능한 라이브러리 확인 
```bash
polyply -list-lib

polyply -list-blocks martini3
```

`gen_params` `gen_itp` `gen_coords` `gen_seq`

- linear polymer 생성
```bash
polyply gen_itp -lib <library_name> -name <name> -seq <monomer:number> -o <name_outfile.itp>
```

- initial structure 생성 #1 
```bash
polyply gen_coords -p <top> -o <name_outfile.gro> -name <name of molecule> -dens <density>
```

- initial structure 생성 #2 
```bash
polyply gen_coords -p <top> -o <name_outfile.gro> -name <name of molecule> -box <x> <y> <z>
```

- 이미 있는 상태에서 구조 생성
```bash
polyply gen_coords -p <top> -o <name_outfile.gro> -name <name of molecule> -c <init_coords.gro> -dens <density>
```

```bash
polyply gen_coords -p <top> -o <name_outfile.gro> -name <name of molecule> -c <init_coords.gro> -box <x> <y> <z>
```


# 2. Tutorial 
### A. GROMOS polymer melts
#### 1) Generating a residue graph
##### A) gen_seq
- `-name` : name of the final molecule
- `-f` : force field files (.ff | .itp)
- `-o` : outfile (.json)
- `-from_string`
	- `<tag>:<number of blocks>:<number of branches>:<resname-probability>` 
	- linear PEG of length 10 : `A:10:1:PEG-1.0`
	- random copolymer of PS-PEG: `A:10:1:PS-0.5,PEG-0.5`
	- branched polymer with 3 generations : `A:3:3:NR3-1`
- `-from_file` 
	- user a molecule defined in an input file as fragment. 
	- format is `<tag:molecule name>`
- `-seq`  
	- Define the sequence order in which to combine macro
	- format is `<macro_tag, macro_tag ...>`
- `-connects` 
	- provide connect records for sequence
	- format is `<seq-index:seq-index:resid-resid>`
	- 만약 첫 번째 block과 3번째 block의 두번째와 세번째 residue로 연결 `<1:3:2-3>`
- `-modf_ter`
	- terminii의 resname을 변경
- `-label` 
	- syntax is `<seqid:label:value-probability,value-probability>`
	- `<0:chiral:R-0.5,S-0.5>`

- 사용법 
```bash
polyply gen_seq -from_string <name:monomers:1:resname-1.0> -o polymer.json -name polymer -seq name -label 0:chiral:R-0.5,S-0.5
```

- atatic PMA 100 repeat unit 
```bash
polyply gen_seq -from_string A:100:1:PMA-1.0 -o PMA.json -name PMA -seq A -label 0:chiral:R-0.5,S-0.5
```

- atatic PMMA 100 repeat unit
```bash
polyply gen_seq -from_string A:100:1:PMMA-1.0 -o PMMA.json -name PMMA -seq A -label 0:chiral:R-0.5,S-0.5
```

- atatic PMMA, PMA 100 repeat unit 
```bash
polyply gen_seq -from_string A:100:1:PMA-0.5,PMMA-0.5 -o PMMA_co_PMA.json -name PMMA_co_PMA -seq A -label 0:chiral:R-0.5,S-0.5
```

#### 2) Generating an itp file
```bash
polyply -list-block 2016H66
```
- [2016H66 논문](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00187)

##### A) gen_params 
  - `-lib` :  force-fields to include from library
  - `-f` :   Input file (ITP|FF)
  - `-o` :  Output ITP (ITP)
  - `-seq` : A linear sequence of residue names.
  - `-seqf` : A graph input file (JSON|TXT|FASTA|IG)

```bash
polyply gen_params -lib <library_name> -o <filename>.itp -name <name of polymer> -seqf <graph-file>.json
```

PMA.ity
```bash
polyply gen_params -lib 2016H66 -o PMA.itp -name PMA -seqf PMA.json
```

PMMA.itp
```bash
polyply gen_params -lib 2016H66 -o PMMA.itp -name PMMA -seqf PMMA.json
```

PMA_co_PMMA.itp
```bash
polyply gen_params -lib 2016H66 -o PMA_co_PMMA.itp -name PMA_co_PMMA -seqf PMMA_co_PMA.json
```

system.top
```
#include "2016H66.ff/forcefield.itp"
#include "PMA.itp"
#include "PMMA.itp"
#include "PMA_co_PMMA.itp"
[ system ]
; name
Melt of PMA, PMMA, PMMA-co-PMA

[ molecules ]
; name  number
PMA 20
PMMA 20
PMA_co_PMMA 10
```

#### 3) Generating starting coordinates

##### A) gen_coords 
- `-p` : topology file (.top)
- `-o` : output gro (.gro)
- `-c` : input file molecules (.gro)
- `-mc` : input file meta-molecule
- `-b` : input file; specify molecule specific building options 
- `-lib` : molecule templates can be loaded from FF 
- `-res` : list of residues to built
- `-ign` : molecules to ignore when building structure
- `-cycles` : molecule names of cyclic molecules
- `-cycle_tol` : tolerance in nm for cycles to be closed
- `-split` : split single residue into more residues
	- format is `<resname><new_resname>-<atom1>,<atom2><etc>:<new_resname>`
- `-lig` : specify ligands of molecules
- `-gs` : grid spacing
- `-grid` : file with grid-points
- `-mi` : max number of tries to grow a molecule
- `-start` : specify which residue to build first
- `-dens` : density of system (kg/m3)
- `-box` : rectangular box x, y, z (nm)

```bash
polyply gen_coords -p <topfilename> -o <filename>.gro -name <name> -dens <density>
```

```bash
polyply gen_coords -p system.top -o melt.gro -name melt -dens 1100
```
- vmd로 melt 구조 확인 

```bash
gmx grompp -f min.mdp -p system.top -c melt.gro -o em.tpr -maxwarn 1

gmx mdrun -v -deffnm em

vmd em.gro

gmx grompp -f pre_NpT.mdp -p system.top -c em.gro -o eq.tpr -maxwarn -1

gmx mdrun -v -deffnm eq
```

### B. Martini Polymer

#### 1). PEO-200 melt 

- 라이브러리 확인
```bash
polyply -list-lib

polyply -list-blocks martini3
```

- gen_seq
```bash
polyply gen_seq -from_string A:200:1:PEO-1.0 -o PEO.json -name PEO -seq A 
```
- gen_params
```bash
polyply gen_params -lib martini3 -o PEO200.itp -name PEO -seqf PEO.json
```

topol.top
```
#include "martini_v3.0.0.itp"
#include "PEO200.itp"
[ system ]
; name
Melt of PEO 200

[ molecules ]
; name  number
PEO 200
```

- gen_coords
```bash
polyply gen_coords -p system.top -o melt.gro -name PEO -dens 1000
```

- minimization & equilibration 
```bash
gmx grompp -f min.mdp -c melt.gro -p topol.top -o min.tpr -maxwarn 1
gmx mdrun -v -deffnm min

gmx grompp -f pre_NpT.mdp -c min.gro -p topol.top -o eq.tpr
gmx mdrun -v -deffnm eq
```

#### 2) PEO_PS blend





#### 3) PS in solvent 




#### 4) PEO_PS copolymer (find problem)

- gen_seq
```bash
polyply gen_seq -from_string A:5:1:PEO-1.0 B:5:1:PS-1.0 -o PEO_PS_PEO.json -name PEO_PS -seq A B A A -connects 0:1:4-0 1:2:4-0 2:3:4-0
```
- gen_seq
```bash
polyply gen_seq -from_string A:5:1:PEO-1.0 B:5:1:PS-1.0 C:10:1:PEO-1.0 -o PEO_PS_PEO.json -name EO_S -seq A B C -connects 0:1:4-0 1:2:4-0
```

- gen_params
```bash
polyply gen_params -lib martini3 -o PEO_PS.itp -name EO_S -seqf PEO_PS_PEO.json
```

topol.top
```
#include "martini_v3.0.0.itp"
#include "PEO_PS.itp"
[ system ]
; name
Melt of EO_S 50

[ molecules ]
; name  number
EO_S 50
```

- gen_coords
```bash
polyply gen_coords -p topol.top -o melt.gro -name EO_S -dens 1000
```

- minimization & equilibration 
```bash
gmx grompp -f min.mdp -c melt.gro -p topol.top -o min.tpr -maxwarn 1
gmx mdrun -v -deffnm min

gmx grompp -f pre_NpT.mdp -c min.gro -p topol.top -o eq.tpr
gmx mdrun -v -deffnm eq
```

- connection 에서 심한 overlap 발생 
```
gmx make_ndx -f eq.tpr 
>>> l

gmx trjconv -f eq.xtc -s eq.tpr -pbc mol -conect -sep -o test_.pdb -n index.ndx
```
- 결과 확인 


### C. PEGylated Lipid Bilayer

#### 1) Generating an itp files 

- 현재 PEGylated Lipid는 DPPE, POPE, DOPE 만 지원
```bash
polyply gen_params -f martini_v3.0.0_phospholipids_v1.itp -lib martini3 -o PEL.itp -name PEL -seq POPE:1 PEO:50 OHter:1
```
- DPPC와 DPPE로 각각 생성해서 비교

#### 2) Generating starting coordinates

##### A) mushroom state
- Generating a bilayer system (95 POPC, 5 POPE )
```bash
insane -l POPC:95 -l POPE:5 -o bilayer.gro -dz 15.0 -x 8.0 -y 8.0 
```

topol.top
```
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "PEL.itp"
[ system ]
; name
PEGylated lipid mushroom

[ molecules ]
; name  number
PEL          5
POPC         95
PEL          5
POPC         95
```


```bash
polyply gen_coords -p topol.top -o struct.gro -res PEO OHter -c bilayer.gro -name test -box 8.0 8.0 20.0 -split POPE:HEAD-NH3:TAIL-PO4,GL1,GL2,C1A,D2A,C3A,C4A,C1B,C2B,C3B,C4B
```
- `-split` : split single residue into more residues
	- format is `<resname><new_resname>-<atom1>,<atom2><etc>:<new_resname>`

- `ValueError: operands could not be broadcast together with shapes (3,) (9,)`
	- bilayer의 box dimension 수정

```bash
polyply gen_coords -p system.top -o struct.gro -res PEO OHter -c bilayer.gro -name test -box 8.0 8.0 20.0 -split POPE:HEAD-NH3:TAIL-PO4,GL1,GL2,C1A,D2A,C3A,C4A,C1B,C2B,C3B,C4B
```

##### B) brush regime

```bash
insane -l POPE:25 -l POPC:75 -o bilayer_brush.gro -dz 33 -x 8.0 -y 8.0 -pbc rectangular
```

topol_brush.top
```
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "PEL.itp"
[ system ]
; name
PEGylated lipid mushroom

[ molecules ]
; name  number
PEL          25
POPC         75
PEL          25
POPC         75
```

build file
```
[ molecule ]
; molname
PEL  0  25
[ rw_restriction ]
; resname start stop  normal   angle
 PEO    3   52    0 0 1     90.0

[ molecule ]
; molname
PEL  100 125
[ rw_restriction ]
; resname start stop  normal   angle
 PEO    2   51    0 0 -1     90.0
```
- `rw_restriction`
	- `<resname> <starting resid> <last resid> <x> <y> <z> <angle in degree>`
```
[ molecule ]
; molname
W  200 14412
[ rectangle ]
; resname start stop  inside-out  x  y  z    a  b  c
W  0  1   out               4.0 4.0 16.5  4.1 4.1 1.5
```

```
[ molecule ]
; molname
NA  200 14412
[ rectangle ]
; resname start stop  inside-out  x  y  z    a  b  c
NA  0  1   out               4.0 4.0 16.5  4.1 4.1 1.5
```

- build file build option 
```
[ template ]
resname <resname>
[ atoms ] 
<atomname> <atomtype> <x> <y> <z>
[ bonds ]
<atomname> <atomname>
```

- `W 12412` `topol_brush.top` 에 추가 
```
polyply gen_coords -p topol_brush.top -o struct.gro -res PEO OHter -c bilayer_brush.gro -name test -box 8.0 8.0 33.0 -split POPE:HEAD-NH3:TAIL-PO4,GL1,GL2,C1A,D2A,C3A,C4A,C1B,C2B,C3B,C4B -b build_file.bld 
```

```bash

gmx grompp -f min.mdp -p topol_brush.top -c struct.gro -o min.tpr -maxwarn 1

gmx mdrun -v -deffnm min

gmx grompp -f pre_NpT.mdp -p topol_brush.top -c min.gro -o eq.tpr

gmx mdrun -v -deffnm eq

gmx make_ndx -f eq.tpr

gmx trjconv -f eq.xtc -s eq.tpr -pbc mol -o sys.pdb -conect -n index.ndx
```

```bash
gmx trjconv -f eq.xtc -s eq.tpr -pbc mol -o eq_pbc.xtc

vmd min.gro eq_pbc.xtc

>>> source cg_bonds_v5.tcl

>>> cg_bonds -top topol_brush.top
```

# 3. Insane 

[공식 링크](https://github.com/Tsjerk/Insane)

### A. 설치 
```
pip install insane
```

### B. options 
- `-f` : Input gro or pdb file (protein)
- `-o` : out gro file 
- `-p` : topology file 
- `-pbc` : PBC type (hexagonal, rectangular, square, cubic) (default: hexagonal)
- `-d` : Distance between periodic image (nm) (default 0)
- `-dz` : Z distance between periodic image (nm)
- `-x` : x dimension (nm)
- `-y` : y dimension (nm)
- `-z` : z dimension (nm)
- `-box` : box in gro or pdb format (3, 9 float for gro, 6 floats for pdb)
- `-n` : index file
- `-l` : lower leaflet lipid type
- `-u` : upper leaflet lipid type
- `-di` : size for inter-leaflet space
- `-a` : Area per lipid (nm2) (default : 0.6 )
- `-au` : Area per lipid for upper leaflet

##### A. Creating a simple lipid bilayer
```bash
insane -l POPC -d 10 -dz 7 -sol W -o bilayer.gro -p topol.top
```

##### B. Creating a more complex bilayer 
```bash
insane -l DPPC:4 -l DIPC:3 -l CHOL:3 -sol W:90 -sol WF:10 -x 20 -y 20 -z 10 -p topol.top -o bilayer.gro
```