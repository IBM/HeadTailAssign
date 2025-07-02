# Head and Tail assignment 

A software that assign the head and tail of polymers.

## Documentation

### **Installation**

* Download the code file to your desired directory and unzip it.

To create the working environment run:


```
conda env create -f environment.yml
```

GAMESS is also a requirement. This software was developed to use GAMESS version 2020-R2 and 2024-R2

### **Running the script**

#### Input file

* The input file provided should contain all the data in one csv file.

* If the user wants to provide the reaction smiles: The csv file should have two columns: 'name', 'reaction'. The 'name' indicates the name of the polymer and the 'reaction' should represent the polymerization reaction. One example of input for reactions:
```
  'name','reaction'
  'polyisobutylene succinic anhydride', '[C]1(=[O])[O][C](=[O])[CH]=[CH]1.[C]1([CH3])[C]([CH3])=[CH][CH]=[CH][CH]=1>>[CH3][C]([CH2][CH]1[C](=[O])[O][C](=[O])[CH2]1)=[CH2]'
```

* If the user wants to provide the monomer: The csv file should have two columns: 'name', 'monomers_list'. The 'name' indicates the name of the polymer and the 'monomers_list' should represent the monomers if homopolymer or list of the monomers separated by '.' if copolymer. One example of input for reactions: One example of input for monomers:
```
'name', 'monomers_list'
'Poly(vinyl formate)','O=COC=C'
'Poly(trimethylene succinate)','OC(=O)CCC(=O)O.OCCCO'
```

* The [main.py](main.py) file has all the steps needed to run the code.
  

To run the code, type at your terminal:

```
python main.py [data]
```
**More information about the functions can be found at the HeadTailAssign module**

#### Output file

* The output file is an csv file with the following structure:
  
```
"name","monomers_list","monomers","polymer_id","classes","mechanism","head_tail"
"Nylon 3 - Poly(propiolactam)","C1CNC1=O","['C1CNC1=O']","polymer_1","polyamide","polycondensation","CN([*:1])C([*:2])=O"
"Poly(hexylene succinate)","OC(=O)CCC(=O)O.OCCCCCCO","['OC(=O)CCC(=O)O', 'OCCCCCCO']","polymer_36","polyester","polycondensation","O=C([*:1])CCC(=O)OCCCCCCO[*:2]"
```

---
## Observations

* Be aware that the Find Monomer method works better with FingerprintSimilarity() metric set to AllBitSimilarity, please change the metric at: 'rdkit\DataStructs\__init__.py'

* Add copolymerizations with products separated by '.'

* Validation tests were performed with theory level HF/STO-3G, but other theory levels can be implemented. To implement a new theory level add a new block of code in the method "_generate_gamess_input_file()" located in [gamess_helper.py](HeadTailAssign/gamess_helper.py). You should modify all the lines shown in the following code to generate the .inp file for GAMESS with the desired theory level:

```
        if (runtype == 'scf' and basis == 'sto3g'):
            output_file.write('!   File created by MacMolPlt 7.7.2 \n')
            output_file.write(' $CONTRL SCFTYP=RHF RUNTYP=ENERGY MAXIT=30 MULT=1 $END \n')
            output_file.write(' $SYSTEM TIMLIM=525600 MWORDS=50 MEMDDI=200 $END \n')
            output_file.write(' $BASIS GBASIS=STO NGAUSS=3 $END \n')
            output_file.write(' $SCF DIRSCF=.TRUE. $END')
            output_file.write('\n')
```

---
## Authorship


* Author: **Brenda Ferrari**
* Co-author: **Ronaldo Giro** 
