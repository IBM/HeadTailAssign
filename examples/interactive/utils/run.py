from IPython.display import SVG, display
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage
import pandas as pd
import io
import os
import subprocess
import shutil
from subprocess import run, Popen, PIPE
import csv
import sys

class Run():
    
    def run_HTA(self, inp):

        if isinstance(inp[0].value, str):
            data = {"name": ["polymer"], "monomers_list": [inp[0].value]}
            df = pd.DataFrame(data)
            display(df)
            df.to_csv("_temp_file.csv", sep=",", quoting=csv.QUOTE_ALL, index=False, header=True)
        else:
            input_file = list(inp[0].value.values())[0]
            content = input_file['content']
            content = io.StringIO(content.decode('utf-8'))
            df = pd.read_csv(content)
            display(df)
            df.to_csv("_temp_file.csv", sep=",", quoting=csv.QUOTE_ALL, index=False, header=True)

        path_to_hta = "..\\..\\"
        temp_path = os.getcwd()

        os.chdir(path_to_hta)
        if os.path.isdir(f"output-directory") == True:
            shutil.rmtree(f"output-directory", ignore_errors=True)

        process = subprocess.Popen(["python", "main.py", "examples\\interactive\\_temp_file.csv"], stdout=subprocess.PIPE)
        while True:
            line = process.stdout.readline()
            if not line:
                break
            print(line.strip().decode('utf-8'))
            sys.stdout.flush()
        df_out = pd.read_csv("results.csv")
        shutil.move("results.csv", "..\\results.csv")
        display(df_out)
        mol = [Chem.MolFromSmiles(m) for m in list(df_out.head_tail)]
                    
        display(MolsToGridImage(mols=mol, legends=list(df_out.name), molsPerRow=3, subImgSize=(400, 400)))
        
        os.chdir(temp_path)
        os.remove('_temp_file.csv')

        
