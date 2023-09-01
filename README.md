# Head and Tail assignment 

A software that assign the head and tail of polymers.

## Documentation

### **Installation**

* Download the code file to your desired directory and unzip it.

To create the working environment run:


```
conda env create
```

GAMESS is also a requirement. This software was developed to use GAMESS version 2020-R2-pgiblas

### **Running the script**

* If the user wants to provide the reaction smiles: The csv file should have two columns: 'name', 'reaction'

* If the user wants to provide the monomer: The csv file should have two columns: 'name', 'monomers_list'

* The [main.py](main.py) file has all the steps needed to run the code.

To run the code, type at your terminal:

```
python main.py [data]
```

**More information about the functions can be found at the HeadTailAssign module**

---
## Observations

Be aware that the Find Monomer method works better with FingerprintSimilarity() metric set to AllBitSimilarity, please change the metric at: 'rdkit\DataStructs\__init__.py'

Add copolymerizations with products separated by '.'

---
## Authorship


* Author: **Brenda Ferrari** ([bferrari](https://github.ibm.com/bferrari))