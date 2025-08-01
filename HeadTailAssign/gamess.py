import glob
import sys
import os
from rdkit import Chem
from openbabel import pybel
from rdkit.Chem import Draw
from tqdm import tqdm
from time import sleep

from HeadTailAssign.gamess_helper import GamessHelper


class GAMESS:
    '''Class related to GAMESS simulation.

    Methods
        generate_gamess_input(self, df, runtype: str, basis: str, FF: str, steps: str) = Generates gamess input files. Generates files, returns None.
        run_gamess_for_files(self) = Run GAMESS simulations for inp files. Returns None.
    '''

    def generate_gamess_input(self, df, name_dir: str, runtype: str, basis: str, FF: str, steps: str):
        '''
        Generates gamess input files. Generates files, returns None.

        Arguments:
            df(dataframe) = The input dataframe generated after mechanisms definition.
            name_dir(str) = Name of the output directory
            runtype(str) = GAMESS input.
            basis(str) = GAMESS input. 
            FF(str) = GAMESS input
            steps = GAMESS input.
        '''
        helper = GamessHelper()
        runtype = runtype  # runtype = scf or opt
        basis = basis

        os.chdir(name_dir)
        
        for name, monomer in tqdm(zip(df['polymer_id'], df['monomers']), total=df.shape[0]):
            sleep(0.1)
            if os.path.isdir(name) == False:
                print(f'\nFile {name} will be created')
                os.mkdir(name)
                os.chdir(name)
            else:
                print(f'\nWARNING: File {name} will be overwritten')
                os.chdir(name)
            m_count = 0
            if isinstance(monomer, list): 
                while m_count < len(monomer):
                    smiles = monomer[m_count]

                    arq_mol = name + '_monomer' + str(m_count) + '.png'
                    arq_molindex = name + '_monomer' + str(m_count) + '-index.png' 

                    mol = Chem.MolFromSmiles(smiles)
                    Draw.MolToFile(mol,arq_mol,size=(300,300))

                    helper._generate_mol_with_atom_index(mol)
                    Draw.MolToFile(mol,arq_molindex,size=(300,300))

                    FF = FF
                    steps = steps
                    molfile = './' + name + '_monomer' + str(m_count) + '.xyz'
                    filename = name + '_monomer' + str(m_count)

                    pbmol = pybel.readstring('smi', smiles)
                    pbmol.make3D(FF)  # make 3D coordinates roughly
                    pbmol.addh()
                    pbmol.localopt(FF, steps)
                    pbmol.write(format='xyz', filename=molfile, overwrite=True)

                    helper._generate_gamess_input_file(name_dir, runtype, basis, filename)

                    m_count+=1
            else:
                smiles = monomer

                arq_mol = name + '_monomer' + str(m_count) + '.png'
                arq_molindex = name + '_monomer' + str(m_count) + '-index.png' 

                mol = Chem.MolFromSmiles(smiles)
                Draw.MolToFile(mol,arq_mol,size=(300,300))

                helper._generate_mol_with_atom_index(mol)
                Draw.MolToFile(mol,arq_molindex,size=(300,300))

                FF = FF
                steps = steps
                molfile = './' + name + '_monomer' + str(m_count) + '.xyz'
                filename = name + '_monomer' + str(m_count)

                pbmol = pybel.readstring('smi', smiles)
                pbmol.make3D(FF)  # make 3D coordinates roughly
                pbmol.addh()
                pbmol.localopt(FF, steps)
                pbmol.write(format='xyz', filename=molfile, overwrite=True)

                helper._generate_gamess_input_file(name_dir, runtype, basis, filename)


            os.chdir("../")

        count = 0
        for file in glob.glob("*.inp"):
            # # Test for WARNING
            # f = open(file, 'r+')
            # f.truncate(0)
            if os.stat(file).st_size >= 100:
                print(f"{file} was generated sucessfully.")
            else:
                count += 1
                print(f"WARNING: {file} could not be generated.")

        if count > 2:
            print("ERROR 04: Too many GAMESS input files could not be generated. Please verify.") # add function responsible
            sys.exit(4)
        
        os.chdir("../")

        return

    def run_gamess_for_files(self, name_dir: str):
        '''
        Run GAMESS simulations for inp files. Returns None.

        Arguments:
            name_dir(str) = Name of the output directory
        '''
        for root, dirs, files in os.walk(name_dir):
            for file in files:
                path = root + '\\' + file # use / to run in linux
                if file.endswith(".inp"):
                    os.system('python ././run_gamess_job.py '+path+' 2')

        count = 0
        for file in glob.glob("*.log"):
            # # Test for WARNING
            # f = open(file, 'r+')
            # f.truncate(0)
            if os.stat(file).st_size >= 100:
                print(f"{file} was generated sucessfully.")
            else:
                count += 1
                print(f"WARNING: {file} could not be generated. Please verify.")

        if count > 2:
            print("ERROR 05: Too many GAMESS output files could not be generated. Please verify.") # add function responsible
            sys.exit(5)

        #os.chdir("../")
        return
