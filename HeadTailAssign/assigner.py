
import struct
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import re
import pandas as pd
import os
import numpy as np
from functools import reduce
from tqdm import tqdm
from time import sleep

from HeadTailAssign.assigner_helper import AssignerHelper


class Assigner:
    '''Class related to assigning head and tail on polymers.

    Methods
        clean_atom_mapping(self, df) -> list = Clean the atom mappings from the reaction smiles. Returns a list.
        create_name(self, df) = Format molecule names to not disturb GAMESS functioning. Returns a dataframe.
        find_monomer(self, df) -> list = Finds the monomer that will be polymerized using chemical similarity. Returns a list.
        get_class(self, df, name_dir) -> list = Get classes of polymers. Returns a list.
        get_mechanism(self, df, name_dir) -> list = Get mechanism of polymers. Returns a list.
        get_head_tail(self, df, name_dir) -> list = Get head and tail of polymers. Returns a list.
    '''
    def __init__(self) -> None:
        '''Initialize the instance of the class Assigner.'''

        self.org_functions = {
            'alkene': '[C:1]=[C:2]',
            'alkyne': '[C:1]#[C:2]',
            'benzene': 'c1ccccc1',
            'primary_amine': '[c,CX4][NH2]',
            'secondary_amine': '[c,CX4][NH]([c,CX4])',
            'tertiary_amine': '[NX3]([c,CX4])([c,CX4])[c,CX4]',
            'aliphatic_alcohol': 'C([OH])',
            'aromatic_alcohol': 'c:c([OH])',
            'eter': 'C[O]C',
            'aliphatic_alkyl_halide': 'C[F,Cl,Br,I]',
            'aromatic_alkyl_halide': 'c:cC([F,Cl,Br,I])=O', #'c:c[F,Cl,Br,I]',
            'thiol': 'C[SH]',
            'carbonyl': 'C(=O)',
            'aldehyde': 'C(C([H])=O)',
            'ketone': '[#6][CX3](=O)[#6]',
            'ester': '[C,c](C([O][C])=O)',
            'carboxilic_acid': '[C,c](C([O])=O)',
            'primary_amide': 'C(C([NH2])=O)',
            'secondary_amide': 'C(C([NH]([CX4]))=O)',
            'tertiary_amide': 'C(C([N]([CX4])([CX4]))=O)',
            'nitrile': 'C([CX2]#[NX1])',
            'dissulfide': 'C[S][S]C',
            'acid_halide': 'C(C([F,Br,Cl,I])=O)',
            'acid_anhydride': 'C(C([O](C(C)=O))=O)',
            'cyanate': '[N]=[C]=[O]',
            'O-heterocycle': 'C@O@C',
            'NC=O-heterocycle': 'C@C([NH](@C))=O'}

    def clean_atom_mapping(self, df) -> list:
        '''
        Clean the atom mappings from the reaction smiles. Returns a list.

        Arguments:
            df(dataframe) = The input dataframe with first column being name (name) and second reaction smiles (reaction).
        '''
        clean_smiles = []
        for row in df.reaction:
            cleanSmiles = re.sub('\:\d{1,3}|\|.+|^\[\'| |\'\]$', '', str(row))
            clean_smiles.append(cleanSmiles)

        return clean_smiles
    
    def create_name(self, df):
        '''
        Format molecule names to not disturb GAMESS functioning. Returns a dataframe.

        Arguments:
            df(dataframe) = The input dataframe with first column being name (name) and second reaction smiles (reaction).
        '''
        df['polymer_id'] = [f'polymer_{i+1}' for i in range(len(df.index))]

        return df 

    def find_monomer(self, df) -> list:
        '''
        Finds the monomer that will be polymerized using chemical similarity. Returns a list.

        Arguments:
            df(dataframe) = The input dataframe with first column being name (name) and second reaction smiles (reaction).
        '''
        maxReactant = [] 
        helper = AssignerHelper()
        name = df.name
        reactantsInp = [helper._separate_reactants(row) for row in df.reaction_clean]
        reactants = reduce(lambda x, y: x+y, reactantsInp)

        productsInp = [helper._separate_products(row) for row in df.reaction_clean]
        products = reduce(lambda x, y: x+y, productsInp)
        for n, x, y in tqdm(zip(name, products, reactants), total=df.shape[0]):
            sleep(0.1)
            productList = [Chem.MolFromSmiles(products) for products in x]
            productFps = [Chem.RDKFingerprint(pl, maxPath=7) for pl in productList]

            dictReact = {}
            reactantList = [Chem.MolFromSmiles(reactants) for reactants in y]
            reactantFps = [Chem.RDKFingerprint(rl, maxPath=7) for rl in reactantList]
            
            if len(productFps) == 1:
                i = 0
                max_sim = 0
                test = []
                while i < len(reactantFps):
                    dictReact[reactantFps[i]] = y[i]

                    compare = DataStructs.FingerprintSimilarity(productFps[0],reactantFps[i])
                    if compare >= 0.0:
                        if compare == 1.0:
                            t = dictReact.get(reactantFps[i])
                            test.append(t)
                            break
                        if compare > max_sim:
                            max_sim = compare
                            if len(test) == 0:
                                t = dictReact.get(reactantFps[i])
                                test.append(t)
                            elif len(test) == 1:
                                del test[0]
                                t = dictReact.get(reactantFps[i])
                                test.append(t)
                    if reactantFps[-1] and len(test) == 0:
                        p = reduce(lambda x, y: x+y, x)
                        test.append(p)
                        break

                    i+=1
                maxReactant.append(test)
            
            if len(productFps) > 1:
                j = 0
                test = []
                while j < len(productFps):
                    i = 0
                    max_sim = 0
                    while i < len(reactantFps):
                        dictReact[reactantFps[i]] = y[i]
                        
                        compare = DataStructs.FingerprintSimilarity(productFps[j],reactantFps[i])
                        if compare >= 0.0:
                            if compare == 1.0:
                                t = dictReact.get(reactantFps[i])
                                test.append(t)
                                break
                            if compare > max_sim:
                                max_sim = compare
                                if len(test) == 0:
                                    t = dictReact.get(reactantFps[i])
                                    test.append(t)
                                elif len(test) == 1:
                                    #del test[0]
                                    t = dictReact.get(reactantFps[i])
                                    test.append(t)
                                elif len(test) == 2:
                                    del test[1]
                                    t = dictReact.get(reactantFps[i])
                                    test.append(t)
                        elif reactantFps[-1] and len(test) == 0:
                            p = reduce(lambda x, y: x+y, x)
                            test.append(p)
                            
                        i+=1
                    j+=1
                maxReactant.append(test)

        allMonomers = []
        for i in maxReactant:
            monomers = list(set(i))
            allMonomers.append(monomers)

        # ERROR 01 test
        # maxReactant.clear()
        if len(allMonomers) == df.shape[0]:
            print("All monomers were extracted.")
        else:
            print("WARNING: Monomers were not generated. Please verify the ouput later.") #add monomer not recognized here
            sys.exit(1)

        return allMonomers

    def get_reactants(self, df) -> list:
        '''
        Non-public method. It gets the reactants from the reaction string.

        Arguments:
            r(str) = Each row from the dataframe which will be iterated.
        '''
        helper = AssignerHelper()
        reactants = [helper._get_reactants(row) for row in tqdm(df.monomers_list, total=df.shape[0])]
        
        return reactants

    def get_class(self, df, name_dir) -> list:
        '''Get classes of polymers. Returns a list.
        
        Arguments:
            df(dataframe): the input dataframe already modified.
            name_dir(str): name of the working directory.
        '''
        os.chdir(name_dir)
        helper = AssignerHelper()
        classes = []

        print("Assigning classes...")
        for name, monomer in tqdm(zip(df['polymer_id'], df['monomers']), total=df.shape[0]):
            sleep(0.1)
            
            if isinstance(monomer, list): 
                if len(monomer) == 1:
                    mol = Chem.MolFromSmiles(monomer[0])

                    org_functions = []
                    atom_maps = []
                    for key, value in self.org_functions.items():
                        match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                        if bool(match_smiles) == True:
                            org_functions.append(key)
                            atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                            atom_maps.append(atom_map)                
                    a = helper._get_classes_one_class(org_functions)
                    if a:
                        classes.append(a)
                    else:
                        try:    
                            file = pd.read_csv(f"{name}\\{name}_monomer0.csv", sep=',')
                            id_number = file['id'].tolist()
                            atom_mappings = helper._sum_one_to_list_lists(atom_maps)
                            a = helper._get_classes_two_class(id_number, atom_mappings, org_functions)
                            classes.append(a)
                        except(FileNotFoundError):
                            print(f"WARNING: File {name}_monomer0.csv was not found. Please verify GAMESS calculation.")
                            classes.append(["error05: GAMESS could not calculate"])
                
                elif len(monomer) > 1:
                    m_count = 0
                    co_atom_maps = []
                    co_org_functions = []
                    while m_count < len(monomer):

                        mol = Chem.MolFromSmiles(monomer[m_count])

                        org_functions = []
                        atom_maps = []
                        for key, value in self.org_functions.items():
                            match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                            if bool(match_smiles) == True:
                                org_functions.append(key)
                                atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                                atom_maps.append(atom_map)
                        if atom_maps:
                            co_atom_maps.append(atom_maps)
                        if org_functions:
                            co_org_functions.append(org_functions)

                        m_count+=1


                    a = helper._get_classes_one_class(co_org_functions)
                    if a:
                        classes.append(a)
                    else:     
                        co_id_number = []
                        m_count = 0
                        while m_count < len(co_org_functions):
                            try:
                                file = pd.read_csv(f"{name}\\{name}_monomer{m_count}.csv", sep=',')
                                id_number = file['id'].tolist()
                                co_id_number.append(id_number)
                            except(FileNotFoundError):
                                print(f"WARNING: File {name}_monomer{m_count}.csv was not found. Please verify GAMESS calculation.")
                                break
                            m_count+=1
                        
                        if len(co_id_number) == len(co_org_functions): 
                            co_atom_mappings = helper._sum_one_to_list_lists(co_atom_maps) #update function
                            a = helper._get_classes_two_class(co_id_number, co_atom_mappings, co_org_functions)
                            classes.append(a)
                        else:
                            classes.append(["error05: GAMESS could not calculate"])
            else:
                mol = Chem.MolFromSmiles(monomer)

                org_functions = []
                atom_maps = []
                for key, value in self.org_functions.items():
                    match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                    if bool(match_smiles) == True:
                        org_functions.append(key)
                        atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                        atom_maps.append(atom_map)                

                a = helper._get_classes_one_class(org_functions)
                if a:
                    classes.append(a)
                else:
                    try:    
                        file = pd.read_csv(f"{name}\\{name}_monomer0.csv", sep=',')
                        id_number = file['id'].tolist()
                        atom_mappings = helper._sum_one_to_list_lists(atom_maps)
                        a = helper._get_classes_two_class(id_number, atom_mappings, org_functions)
                        classes.append(a)
                    except(FileNotFoundError):
                        print(f"WARNING: File {name}_monomer0.csv was not found. Please verify GAMESS calculation.")
                        classes.append(["error05: GAMESS could not calculate"])

        classes_list = []
        for c in classes:
            if c:
                classes_list.append(c)
            else:
                classes_list.append(['not recognized'])
        classes = helper._reduce_list(classes_list)
        os.chdir('../')
        return classes

    def get_mechanism(self, df, name_dir) -> list:
        '''Get mechanism of polymers. Returns a list.
        
        Arguments:
            df(dataframe): the input dataframe already modified.
            name_dir(str): name of the working directory.
        '''
        if 'reaction_clean' in df:
            os.chdir(name_dir)
            helper = AssignerHelper()
            mechanism = []

            print("Assigning mechanism...")
            for name, reaction_clean, classes in tqdm(zip(df['polymer_id'], df['reaction_clean'], df['classes']), total=df.shape[0]):
                sleep(0.1)

                if classes == 'vinyl':
                    mechanism_vinyl = helper._get_vinyl_mechanism(reaction_clean)
                    mechanism.append(mechanism_vinyl[0])
                elif classes == 'polyamide':
                    mechanism_polyamide = helper._get_polyamide_mechanism(classes)
                    mechanism.append(mechanism_polyamide)
                elif classes == 'polyester':
                    mechanism_polyester = helper._get_polyester_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                elif classes == 'polyurethane':
                    mechanism_polyester = helper._get_polyurethane_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                elif classes == 'polyether':
                    mechanism_polyester = helper._get_polyether_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                else:
                    mechanism.append('not recognized')        

            mechanism_reduce = []
            for m in mechanism:
                if type(m) is list:
                    mechanism_reduce1 = reduce(lambda x, y: x+y, m)
                    mechanism_reduce.append(mechanism_reduce1)
                else:   
                    mechanism_reduce.append(m)
            os.chdir('../')

        else:
            os.chdir(name_dir)
            helper = AssignerHelper()
            mechanism = []
            
            for name, classes in zip(df['polymer_id'], df['classes']):
                
                if classes == 'vinyl':
                    print("WARNING: There is no way to assign the type of vinyl mechanism, please provide the reaction for this functionality")
                    mechanism.append('polyaddition')
                elif classes == 'polyamide':
                    mechanism_polyamide = helper._get_polyamide_mechanism(classes)
                    mechanism.append(mechanism_polyamide)
                elif classes == 'polyester':
                    mechanism_polyester = helper._get_polyester_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                elif classes == 'polyurethane':
                    mechanism_polyester = helper._get_polyurethane_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                elif classes == 'polyether':
                    mechanism_polyester = helper._get_polyether_mechanism(classes)
                    mechanism.append(mechanism_polyester)
                else:
                    mechanism.append(['not recognized'])        

            mechanism_reduce = []
            for m in mechanism:
                if type(m) is list:
                    mechanism_reduce1 = reduce(lambda x, y: x+y, m)
                    mechanism_reduce.append(mechanism_reduce1)
                else:   
                    mechanism_reduce.append(m)
            os.chdir('../')

        return mechanism_reduce 

    def get_head_tail(self, df, name_dir, notation='usual', structure='monomer') -> list:
        '''Get head and tail of polymers. Returns a list.
        
        Arguments:
            df(dataframe): the input dataframe already modified.
            name_dir(str): name of the working directory.
            notation(str): output format. If 'usual' outputs usual notation for ligands (*:1 and *:2)
            structure(str): output format. May be monomer or oligomer.
        '''
        os.chdir(name_dir)
        helper = AssignerHelper()
        head_tail = []

        print("Assigning head and tail...")
        for name, monomer, classes in tqdm(zip(df['polymer_id'], df['monomers'], df['classes']), total=df.shape[0]):
            sleep(0.1)
            if isinstance(monomer, list): 
                if len(monomer) == 1:

                    mol = Chem.MolFromSmiles(monomer[0])

                    [a.SetAtomMapNum(i+1) for i,a in enumerate(mol.GetAtoms())]

                    smiles = (Chem.MolToSmiles(mol))

                    org_functions = []
                    atom_maps = []
                    for key, value in self.org_functions.items():
                        match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                        if bool(match_smiles) == True:
                            org_functions.append(key)
                            atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                            atom_maps.append(atom_map)

                    ##### GET ATOM MAPPING ######
                    smiles_head_tail = helper._get_atom_mapping(smiles, atom_maps, org_functions, classes)
                    if smiles_head_tail:
                        if type(smiles_head_tail) == 'list':
                            smiles_head_tail_reduce = reduce(lambda x, y: x+y, smiles_head_tail)
                            head_tail.append(smiles_head_tail_reduce)
                        else:
                            head_tail.append(smiles_head_tail)
                    else:
                        head_tail.append('not recognized')

                elif len(monomer) > 1:
                    m_count = 0
                    co_head_tail = []
                    co_atom_maps = []
                    co_org_functions = []
                    while m_count < len(monomer):

                        mol = Chem.MolFromSmiles(monomer[m_count])

                        [a.SetAtomMapNum(i+1) for i,a in enumerate(mol.GetAtoms())]

                        smiles = (Chem.MolToSmiles(mol))

                        org_functions = []
                        atom_maps = []
                        for key, value in self.org_functions.items():
                            match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                            if bool(match_smiles) == True:
                                org_functions.append(key)
                                atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                                atom_maps.append(atom_map)
                        co_atom_maps.append(atom_maps)
                        co_org_functions.append(org_functions)

                        ##### GET ATOM MAPPING ######
                        smiles_head_tail = helper._get_atom_mapping(smiles, atom_maps, org_functions, classes)
                        if smiles_head_tail:
                            smiles_head_tail_reduce = reduce(lambda x, y: x+y, smiles_head_tail)
                            co_head_tail.append(smiles_head_tail_reduce)
                        else:
                            co_head_tail.append('not recognized')
                        
                        m_count+=1

                    if co_head_tail:
                        head_tail.append(co_head_tail)
            else:
                mol = Chem.MolFromSmiles(monomer)

                [a.SetAtomMapNum(i+1) for i,a in enumerate(mol.GetAtoms())]

                smiles = (Chem.MolToSmiles(mol))

                org_functions = []
                atom_maps = []
                for key, value in self.org_functions.items():
                    match_smiles = mol.GetSubstructMatches(Chem.MolFromSmarts(value))

                    if bool(match_smiles) == True:
                        org_functions.append(key)
                        atom_map = [element for tupl in match_smiles for element in tupl] #flattening a nested structure.
                        atom_maps.append(atom_map)

                ##### GET ATOM MAPPING ######
                smiles_head_tail = helper._get_atom_mapping(smiles, atom_maps, org_functions, classes)
                if smiles_head_tail:
                    smiles_head_tail_reduce = reduce(lambda x, y: x+y, smiles_head_tail)
                    head_tail.append(smiles_head_tail_reduce)
                else:
                    head_tail.append('not recognized')

        head_tail_list = []
        for h in head_tail:
            if h:
                head_tail_list.append(h)
            else:
                head_tail_list.append('not recognized')
        
        if notation == 'usual':
            head_tail_list = self._change_headtail_notation(head_tail_list)
        if structure == 'oligomer':
            head_tail_list = self._create_oligomer(df, head_tail_list)

        os.chdir('../')
        
        return head_tail_list
    
    def _change_headtail_notation(self, head_tail_list) -> list:
        '''Change head and tail to usual notation. Returns a list.
        
        Arguments:
            head_tail_list(list): the list of monomers with head and tail assigned.
        '''
        new_notation = []

        for monomer in head_tail_list:
            if type(monomer) == list:
                try:
                    new_notation_monomers = []
                    for m in monomer:
                        m_h = re.sub('\:\*h', '*:1', m)
                        m_t = re.sub('\:\*t', '*:2', m_h)
                        new_notation_monomers.append(m_t)
                    
                    new_notation.append(new_notation_monomers)
                except: 
                    new_notation.append(f'WARNING: not translated - {monomer}')
            else:
                try:
                    m_h = re.sub('\:\*h', '*:1', m)
                    m_t = re.sub('\:\*t', '*:2', m_h)
                    new_notation.append(m_t)
                except:
                    new_notation.append(f'WARNING: not translated - {monomer}')

        return new_notation

    def _create_oligomer(self, df, head_tail_list) -> list:
        '''Create oligomers from the list of monomers. Returns a list.
        
        Arguments:
            head_tail_list(list): the list of monomers with head and tail assigned.
        '''
        oligomer_list = []
        for name, monomer, classes in zip(df['polymer_id'], head_tail_list, df['classes']):

            if type(monomer) == list:
                try:
                    monomers = []
                    if classes == 'polyamide':
                        for m in monomer:
                            m_rn_h: str = re.sub('\*:1', 'XeH', m)
                            m_rn_tn = re.sub('\*:2', 'Rn', m_rn_h)
                            if bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('C([Rn])(=O)'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(0, Chem.MolFromSmiles(m_rn_tn))

                            elif bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('N[XeH]'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(1, Chem.MolFromSmiles(m_rn_tn))
                        smarts = AllChem.ReactionFromSmarts('([C:1]([Rn])=[O:2]).[N:3][XeH]>>([C:1]([N:3])=[O:2])')

                    if classes == 'polyester':
                        for m in monomer:
                            m_rn_h: str = re.sub('\*:1', 'XeH', m)
                            m_rn_tn = re.sub('\*:2', 'Rn', m_rn_h)

                            if bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('C([Rn])(=O)'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(0, Chem.MolFromSmiles(m_rn_tn))

                            if bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('O[XeH]'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(1, Chem.MolFromSmiles(m_rn_tn))
                        smarts = AllChem.ReactionFromSmarts('([C:1]([Rn])=[O:2]).[O:3][XeH]>>[C:1](=[O:2])[O:3]')

                    elif classes == 'polyurethane':
                        for m in monomer:
                            m_rn_h: str = re.sub('\*:1', 'XeH', m)
                            m_rn_tn = re.sub('\*:2', 'Rn', m_rn_h)

                            if bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('[N](C([Rn])=O)'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(0, Chem.MolFromSmiles(m_rn_tn))

                            if bool(Chem.MolFromSmiles(m_rn_tn).GetSubstructMatches(Chem.MolFromSmarts('CO[XeH]'))) == True:
                                if monomers == []:
                                    monomers.append(Chem.MolFromSmiles(m_rn_tn))
                                else:
                                    monomers.insert(0, Chem.MolFromSmiles(m_rn_tn))
                        smarts = AllChem.ReactionFromSmarts('[N:1]([C:2]([Rn])=[O:3]).[C:4][O:5][XeH]>>[N:1]([C:2]([O:5][C:4])=[O:3])')

                    oligomer = smarts.RunReactants(monomers)

                    result = []
                    for i in range(len(oligomer)):
                        m_r_h = Chem.MolToSmiles(oligomer[i][0])
                        m_h: str = re.sub('XeH', '*:1', m_r_h)
                        mol_monomer = re.sub('Rn', '*:2', m_h)
                        if mol_monomer not in result:
                            result.append(mol_monomer)
                    oligomer_list.append(result)
                except:
                    oligomer_list.append(monomer)

            else:
                oligomer_list.append(monomer)
                
        oligomer_result = [element for tupl in oligomer_list for element in tupl]

        return oligomer_result