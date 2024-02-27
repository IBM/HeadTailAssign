from itertools import product
from rdkit import Chem
from rdkit import DataStructs
import re
from functools import reduce
from time import sleep


class AssignerHelper:
    '''
    Class related to Non-public methods which are Assigner class helpers.

    Methods
        _separate_reactants(self, row) = Non-public method. It separates the reactants from the reaction string.
        _separate_catalysts(self, row) = Non-public method. It separates the catalysts from the reaction string.
        _separate_products(self, row) = Non-public method. It separates the products from the reaction string.
        _get_reactants(self, r) = Non-public method. It gets the reactants from the reaction string.
        _sum_one_to_list_lists(self, atom_maps)= Non-public method. It adds a value of one to a list of lists.
        _reduce_list(self, classes_list)= Non-public method. It reduces a list of lists depending on object type.

        _get_vinyl_mechanism(self, r, c) = Non-public method. It gets the type of vinyl mechanism involved in the reaction.
        _get_polyamide_mechanism(self, c) = Non-public method. It gets the polyamide mechanism involved in the reaction.
        _get_polyester_mechanism(self, c)= Non-public method. It gets the polyester mechanism involved in the reaction.
        _get_polyurethane_mechanism(self, c) = Non-public method. It gets the polyurethane mechanism involved in the reaction.
        _get_polyether_mechanism(self, c) = Non-public method. It gets the polyether mechanism involved in the reaction.

        _get_classes_one_class(self, org_function)= Non-public method. Gets the classification of the polymerization mechanism for one functional group.
        _get_classes_two_class(self, id_number, atom_map, org_functions) = Non-public method. Gets the classification of the polymerization mechanism for more than one functional group.
        _open_ring(self, classes, smiles) = Non-public method. Open ring of cyclic monomers.
        _sanitize_polymer(self, classes, smiles) = Non-public method. Change the smiles format after opening ring.
        _get_atom_mapping(self, smiles, atom_mapping, org_functions, classes)= Non-public method. Gets the classification of the polymerization mechanism for more than one functional group.
    '''

    def _canonicalize(self, monomers):
        '''
        Non-public method. It canonicalizes the SMILES.

        Arguments:
            monomers(list) = List of monomers.
        '''
        monomers_can_list = []
        if monomers == list():
            for m in monomers:
                monomers_can = Chem.MolToSmiles(Chem.MolFromSmiles(m))
                monomers_can_list.append(monomers_can)
        else:
            monomers_can = Chem.MolToSmiles(Chem.MolFromSmiles(monomers[0]))
            monomers_can_list.append(monomers_can)           

        return monomers_can_list

    def _separate_reactants(self, row):
        '''
        Non-public method. It separates the reactants from the reaction string.

        Arguments:
            row(str) = Each row from the dataframe which will be iterated.
        '''
        reactantsList = []
        if bool(re.findall('>>', str(row))) == True:
            reaction = row.split('>>')

            reactants = reaction[0].split('.')
            reactantsList.append(reactants)
        else:
            reaction = re.findall('.*?(?=>)', str(row))

            reactants = reaction[0].split('.')
            reactantsList.append(reactants)

        return reactantsList

    def _separate_catalysts(self, row):
        '''
        Non-public method. It separates the catalysts from the reaction string.

        Arguments:
            row(str) = Each row from the dataframe which will be iterated.
        '''
        catalystsList = []
        if bool(re.findall('>>', str(row))) == False:
            reaction = re.findall('(?<=>).*?(?=>)', str(row))

            catalysts = reaction[0].split('.')
            catalystsList.append(catalysts)            
            
        return catalystsList

    def _separate_products(self, row):
        '''
        Non-public method. It separates the products from the reaction string.

        Arguments:
            row(str) = Each row from the dataframe which will be iterated.
        '''
        productList = []
        if bool(re.findall('>>', str(row))) == True:
            reaction = row.split('>>')
            if bool(re.findall('.', str(reaction[1]))) == True:
                products = reaction[1].split('.')
                productList.append(products)
            else:
                productList.append(reaction[1])

        else:
            reaction = row.split('>')
            if bool(re.findall('.', str(reaction[2]))) == True:
                products = reaction[2].split('.')
                productList.append(products)
            else:
                productList.append(reaction[2])

        return productList

    def _get_reactants(self, r):
        '''
        Non-public method. It gets the reactants from the reaction string.

        Arguments:
            r(str) = Each row from the dataframe which will be iterated.
        '''
        sleep(0.1)
        reactants = r.split('.')
        
        return reactants 

    def _sum_one_to_list_lists(self, atom_maps):
        '''
        Non-public method. It adds a value of one to a list of lists.

        Arguments:
            atom_maps(list): atom mappings obtained on the get_class method.
        '''        
        if bool(re.findall('^\[\[\[', str(atom_maps))) == True:
            atom_mappingss = []
            for aaa in atom_maps:
                atom_mappings = []
                atom_mappingss.append(atom_mappings)
                for aa in aaa:
                    atom_mapping = []
                    atom_mappings.append(atom_mapping)
                    for a in aa:
                        aaa = a + 1
                        atom_mapping.append(aaa)
        else:
            atom_mappingss = []
            for aa in atom_maps:
                atom_mappings = []
                atom_mappingss.append(atom_mappings)
                for a in aa:
                    aaa = a + 1
                    atom_mappings.append(aaa)

        return atom_mappingss

    def _reduce_list(self, classes_list):
        '''
        Non-public method. It reduces a list of lists depending on object type.

        Arguments:
            classes_list(list): classes obtained by get_class method.
        '''
        classes_reduce_list = reduce(lambda x, y: x+y, classes_list)
        classes = []
        for c in classes_reduce_list:
            if type(c) is list:
                classes_reduce = reduce(lambda x, y: x+y, c)
                classes.append(classes_reduce)
            else:   
                classes.append(c)
        
        return classes

    def _get_vinyl_mechanism(self, row):
        '''
        Non-public method. It gets the type of vinyl mechanism involved in the reaction.

        Arguments:
            r(list) = reactant list
            c(str) = classes that comes from iteration in dataframe.
        '''
        #https://www.tcichemicals.com/assets/brochure-pdfs/Brochure_F2037_E.pdf
        smiles = {
            'anionic': ['CCCC[Li]', '[H]C(=O)c1ccccc1C(=O)N2CCCCC2', 'CC(C(=O)O[Na])c3ccc2oc1ccccc1c(=O)c2c3',
                        'C2CN=C1NCCCN1C2', 'C2CN=C1CCCN1C2', 'C2CCC1=NCCCN1CC2',
                        'CC(=NOC(=O)c1ccccc1)c2ccccc2', 'O=C(NC1CCCCC1)OCc2ccccc2N(=O)=O',
                        'COc3ccc(C(=O)C(OC(=O)NC1CCCCC1)c2ccc(OC)cc2)cc3', 'CC(=O)c1ccccc1C(=O)N2CCCCC2',
                        'COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c2ccccc2N(=O)=O'],
            'cationic': ['FB(F)F', 'c2ccc([I+]c1ccccc1)cc2','O=S(=O)([O-])C(F)(F)F',
                         'F[P-](F)(F)(F)(F)F', 'F[As](F)(F)(F)(F)F', 'CC(C)(C)c2ccc([I+]c1ccc(C(C)(C)C)cc1)cc2',
                         'CC(C)c2ccc([I+]c1ccc(C(C)C)cc1)cc2',
                         'Fc4c(F)c(F)c([B-](c1c(F)c(F)c(F)c(F)c1F)(c2c(F)c(F)c(F)c(F)c2F)c3c(F)c(F)c(F)c(F)c3F)c(F)c4F',
                         'c3ccc([S+](c1ccccc1)c2ccccc2)cc3', '[Br-]', 'F[B-](F)(F)F', 
                         'O=S(=O)([O-])C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'Cc3ccc([S+](c1ccc(C)cc1)c2ccc(C)cc2)cc3',
                         'N#[N+]c1ccc(N(=O)=O)cc1', 'COc2ccc(c1nc(C(Cl)(Cl)Cl)nc(C(Cl)(Cl)Cl)n1)cc2',
                         'Cc2ccc(C=Cc1nc(C(Cl)(Cl)Cl)nc(C(Cl)(Cl)Cl)n1)o2', 'ClC(Cl)(Cl)c2nc(C=Cc1ccco1)nc(C(Cl)(Cl)Cl)n2',
                         'COc2ccc(C=Cc1nc(C(Cl)(Cl)Cl)nc(C(Cl)(Cl)Cl)n1)cc2OC', 'COc2ccc(C=Cc1nc(C(Cl)(Cl)Cl)nc(C(Cl)(Cl)Cl)n1)cc2',
                         'ClC(Cl)(Cl)c3nc(c2ccc1OCOc1c2)nc(C(Cl)(Cl)Cl)n3'],
            'radical': ['CC(C)(C)OO', 'CC1(C)C2CCC1(C)C(=O)C2=O', 'CC(=O)c1ccccc1',
                        'Cc1ccccc1C(=O)c2ccccc2', 'Cc2cccc(C(=O)c1ccccc1)c2',
                        'O=C(c1ccccc1)c2ccccc2', 'CC(=O)c1cccc(O)c1', 'CC(=O)c1ccc(O)cc1',
                        'Cc2ccc(C(=O)c1ccccc1)cc2C', 'O=C(c1ccccc1)c2cccc(O)c2', 
                        'O=C(c1ccccc1)c2ccc(O)cc2', 'O=C(c1ccc(O)cc1)c2ccc(O)cc2',
                        'O=C(O)c2ccc(C(=O)c1ccccc1)cc2', 'CN(C)c2ccc(C(=O)c1ccc(N(C)C)cc1)cc2',
                        'CN(C)c2ccc(C(=O)c1ccccc1)cc2', 'O=C(c2ccc1c(=O)oc(=O)c1c2)c4ccc3c(=O)oc(=O)c3c4',
                        'COC(=O)c1ccccc1C(=O)c2ccccc2', 'O=C(O)c1ccccc1C(=O)c2ccccc2',
                        'CCN(CC)c2ccc(C(=O)c1ccc(N(CC)CC)cc1)cc2', 'O=C(c1ccc(Cl)cc1)c2ccc(Cl)cc2',
                        'O=C(c1ccccc1)c3ccc(c2ccccc2)cc3', 'O=C(c1ccccc1)c3ccc(C(=O)c2ccccc2)cc3',
                        'Cc3ccc(Sc2ccc(C(=O)c1ccccc1)cc2)cc3', 'O=c2c(=O)c1ccccc1c3ccccc23',
                        'COC(=O)C(=O)c1ccccc1', 'COc2ccc(C(=O)C(=O)c1ccc(OC)cc1)cc2',
                        'O=C(C(=O)c1ccccc1)c2ccccc2', 'O=c2c1ccccc1ccc3ccccc23', 'CC(C)(O)C(=O)c1ccccc1',
                        'CC(C)(O)C(=O)c1cccc(OCCO)c1', 'O=C(c1ccccc1)C2(O)CCCCC2', 'O=C(c1ccccc1)C(O)c2ccccc2',
                        'COc2ccc(C(=O)C(O)c1ccc(OC)cc1)cc2', 'CCOC(OCC)C(=O)c1ccccc1', 
                        'CCC(CC)COC(C(=O)c1ccccc1)c2ccccc2', 'CCC(CC)OC(C(=O)c1ccccc1)c2ccccc2',
                        'CCOC(C(=O)c1ccccc1)c2ccccc2', 'COC(C(=O)c1ccccc1)c2ccccc2', 'COC(OC)(C(=O)c1ccccc1)c2ccccc2',
                        'CSc2ccc(C(=O)C(C)(C)N1CCOCC1)cc2', 'CCC(Cc1ccccc1)(C(=O)c3ccc(N2CCOCC2)cc3)N(C)C',
                        'CC(=NO)C(=O)c1ccccc1', 'CC(C)c3ccc2sc1ccccc1c(=O)c2c3',
                        'CCCOc2ccc(Cl)c3C(=O)c1ccccc1Cc23', 'O=c2c1ccccc1c(=O)c3cc(S(=O)(=O)O[Na])ccc23.[H]O[H]',
                        'CCc3ccc2c(=O)c1ccccc1c(=O)c2c3', 'CCc3cc(CC)c2sc1ccccc1c(=O)c2c3',
                        'COc3ccc2sc1ccc(OC)cc1c(=O)c2c3', 'Clc1ccccc1c4nc(c2ccccc2)c(c3ccccc3)n4C8(c5ccccc5Cl)N=C(c6ccccc6)C(c7ccccc7)=N8',
                        'Cc3cc(C)c(C(=O)P(=O)(c1ccccc1)c2ccccc2)c(C)c3', 'Cc3cc(C)c(C(=O)P(=O)(C(=O)c1c(C)cc(C)cc1C)c2ccccc2)c(C)c3',
                        'Cc2cc(C)c(C(=O)P(=O)(O[Li])c1ccccc1)c(C)c2', '[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Fe+2]']
        }
        
        mechanism_vinyl = []
 
        reactants = self._separate_reactants(row)
        reactants_reduce = reduce(lambda x, y: x+y, reactants)
        for r in reactants_reduce:
            mol = Chem.MolFromSmiles(r)
            mol_bits = Chem.RDKFingerprint(mol)
            for value in smiles['anionic']:
                value_mol = Chem.MolFromSmiles(value)
                value_bits = Chem.RDKFingerprint(value_mol)
                if DataStructs.FingerprintSimilarity(mol_bits,value_bits) == 1.0:
                    mechanism_vinyl.append('anionic polymerization')

            for value in smiles['cationic']:
                value_mol = Chem.MolFromSmiles(value)
                value_bits = Chem.RDKFingerprint(value_mol)
                if DataStructs.FingerprintSimilarity(mol_bits,value_bits) == 1.0:
                    mechanism_vinyl.append('cationic polymerization') 

            for value in smiles['radical']:
                value_mol = Chem.MolFromSmiles(value)
                value_bits = Chem.RDKFingerprint(value_mol)
                if DataStructs.FingerprintSimilarity(mol_bits,value_bits) == 1.0:
                    mechanism_vinyl.append('radical polymerization')

        catalysts = self._separate_catalysts(row)
        if catalysts:
            catalysts_reduce = reduce(lambda x, y: x+y, catalysts)
            for c in catalysts_reduce:
                a = Chem.MolFromSmiles(c)
                for value in smiles['anionic']:
                    if a.HasSubstructMatch(Chem.MolFromSmarts(value)):
                        mechanism_vinyl.append('anionic polymerization')

                for value in smiles['cationic']:
                    if a.HasSubstructMatch(Chem.MolFromSmarts(value)):
                        mechanism_vinyl.append('cationic polymerization') 

                for value in smiles['radical']:
                    if a.HasSubstructMatch(Chem.MolFromSmarts(value)):
                        mechanism_vinyl.append('radical polymerization') 
                     
        if not mechanism_vinyl:
            mechanism_vinyl.append('polyaddition')

        return mechanism_vinyl

    def _get_polyamide_mechanism(self, c):
        '''
        Non-public method. It gets the polyamide mechanism involved in the reaction.

        Arguments:
            c(str) = classes that comes from iteration in dataframe.
        '''
        mechanism_polyamide = []
        if c == 'polyamide':
            mechanism_polyamide.append('polycondensation')

        return mechanism_polyamide

    def _get_polyester_mechanism(self, c):
        '''
        Non-public method. It gets the polyester mechanism involved in the reaction.

        Arguments:
            c(str) = classes that comes from iteration in dataframe.
        '''
        mechanism_polyester = []
        if c == 'polyester':
            mechanism_polyester.append('polycondensation')

        return mechanism_polyester

    def _get_polyurethane_mechanism(self, c):
        '''
        Non-public method. It gets the polyurethane mechanism involved in the reaction.

        Arguments:
            c(str) = classes that comes from iteration in dataframe.
        '''
        mechanism_polyurethane = []
        if c == 'polyurethane':
            mechanism_polyurethane.append('polycondensation')

        return mechanism_polyurethane

    def _get_polyether_mechanism(self, c):
        '''
        Non-public method. It gets the polyether mechanism involved in the reaction.

        Arguments:
            c(str) = classes that comes from iteration in dataframe.
        '''
        mechanism_polyether = []
        if c == 'polyether':
            mechanism_polyether.append('polycondensation')

        return mechanism_polyether

    def _get_classes_one_class(self, org_function):
        '''
        Non-public method. Gets the classification of the polymerization mechanism for one functional group.
        
        Arguments:
            org_function(list): organic functions obtained on the get_class method.
        '''
        print(org_function)
        classes = []
        if bool(re.findall('^\[\[', str(org_function))) == False:
            if org_function == ['alkene']:
                classes.append('vinyl')
            elif org_function == ['alkyne']:
                classes.append('vinyl')
            elif org_function == ['primary_amine', 'aliphatic_alkyl_halide', 'acid_halide']:
                classes.append('polyamide')
            elif org_function == ['primary_amine', 'aliphatic_alcohol', 'carboxilic_acid']:
                classes.append('polyamide')
            elif org_function == ['secondary_amide', 'NC=O-heterocycle']: 
                classes.append('polyamide')
            elif org_function == ['ester', 'aliphatic_alcohol', 'carboxilic_acid']:
                classes.append('polyester') 
            elif org_function == ['aliphatic_alcohol', 'carboxilic_acid']:
                classes.append('polyester') 
            elif org_function == ['aliphatic_alcohol', 'cyanate']:
                classes.append('polyurethane') 
            elif org_function == ['eter', 'O-heterocycle']: 
                classes.append('polyether')
            elif org_function == ['carbonyl']:
                classes.append('polyether')
        # EDIT HERE FOR OTHER CLASSES
        else:
            if org_function == [['alkene'], ['alkene']]:
                classes.append('vinyl')
            if org_function == [['primary_amine'], ['aliphatic_alcohol', 'carboxilic_acid']]:
                classes.append('polyamide')
            if org_function == [['aliphatic_alcohol'], ['aliphatic xc_alcohol', 'carboxilic_acid']]:
                classes.append('polyester')

        return classes

    def _get_classes_two_class(self, id_number, atom_map, org_functions):
        '''
        Non-public method. Gets the classification of the polymerization mechanism for more than one functional group.
        
        Arguments:
            id_number(list): id numbers obtained by GAMESS calculation located on csv files.
            org_functions(list): organic functions obtained on the get_class method.
            atom_map(list): atom mappings obtained on the get_class method.
        '''
        # EDIT HERE FOR OTHER CLASSES
        co_id_number = []
        co_id_number.append(id_number[1])
        co_id_number.append(id_number[0])
        if bool(re.findall('^\[\[', str(org_functions))) == True:

            classes_dict = {'vinyl':[[['alkene'], ['alkene']], [['alkene'], ['alkene']]],#
                    'vinyl_A':[[['alkyne'], ['alkyne']], [['alkyne'], ['alkyne']]],#
                    'polyamide':[[['primary_amine'], ['aliphatic_alcohol', 'carboxilic_acid']]], 
                                #[['aliphatic_alcohol', 'carboxilic_acid'], ['primary_amine']]],#
                    'polyamide_A':[[['primary_amine'], ['aliphatic_alkyl_halide', 'acid_halide']], 
                                #[['aliphatic_alkyl_halide', 'acid_halide'], ['primary_amine']], 
                                [['primary_amine'], ['aliphatic_alkyl_halide']]],
                    'polyester':[[['aliphatic_alcohol'], ['aliphatic_alcohol', 'carboxilic_acid']], 
                                #[['aliphatic_alcohol', 'carboxilic_acid'], ['aliphatic_alcohol']],
                                [['carboxilic_acid'], ['aliphatic_alcohol', 'carboxilic_acid'], ['aliphatic_alcohol']],
                                [['benzene', 'aliphatic_alcohol', 'carboxilic_acid'], ['aliphatic_alcohol']]],#
                    'polyurethane':[[['aliphatic_alcohol'], ['cyanate']],
                                    [['eter'], ['cyanate'], ['aliphatic_alcohol']]],
            }
            classes_final = []
            classes = []
            
            for i in co_id_number:
                functions = []
                for j in i:
                    print(f'functions {functions}')
                    for atom_mapping, o_f in zip(atom_map, org_functions):
                        functions_n = []
                        functions_e = []
                        for a_m, oo_f in zip(atom_mapping, o_f):                            
                            if j in a_m:
                                if oo_f == 'alkene':
                                    o_f_ = 'alkene'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'alkyne':
                                    o_f_ = 'alkyne'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'primary_amine':
                                    o_f_ = 'primary_amine'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'acid_halide':
                                    o_f_ = 'acid_halide'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'aliphatic_alkyl_halide':
                                    o_f_ = 'aliphatic_alkyl_halide'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'aliphatic_alcohol':
                                    o_f_ = 'aliphatic_alcohol'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'carboxilic_acid':
                                    o_f_ = 'carboxilic_acid'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'secondary_amide':
                                    o_f_ = 'secondary_amide'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)
                                elif oo_f == 'NC=O-heterocycle':
                                    o_f_ = 'NC=O-heterocycle'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)                  
                                elif oo_f == 'cyanate':
                                    o_f_ = 'cyanate'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)  
                                elif oo_f == 'eter':
                                    o_f_ = 'eter'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e) 
                                elif oo_f == 'O-heterocycle':
                                    o_f_ = 'O-heterocycle'
                                    if a_m in atom_map[0]:
                                        functions_n.append(o_f_)
                                    elif a_m not in atom_map[0]:
                                        functions_e.append(o_f_)
                                    functions.append(functions_n)
                                    functions.append(functions_e)

                                if ['alkene'] not in functions:
                                    functions_a = []
                                    for elem in functions:
                                        if elem not in functions_a:
                                            functions_a.append(elem)
                                    k = functions_a
                                    functions = functions_a
                                
                                functions = [x for x in functions if x != []]
                                #print(f'functions {functions}')
                                for key, value in classes_dict.items():
                                    for v in value:

                                        if sorted(v) == sorted(functions):
                                            if key == 'polyamide_A':
                                                classes.append('polyamide')
                                                break 
                                            else:
                                                print("HERE")           
                                                classes.append(key)
                                                break

                        if len(classes) == 1:
                            break
                if len(classes) == 1:
                    classes_final.append(classes)
                elif len(classes) >= 2:
                    classes_final.append(classes[0])

                if len(classes_final) == 1:
                    break
                            
        else:
            classes_dict = {'vinyl':['alkene'],
                    'vinyl_A':['alkyne'],
                    'polyamide':['primary_amine', 'aliphatic_alcohol', 'carboxilic_acid'],
                    'polyamide_A':['primary_amine', 'aliphatic_alkyl_halide', 'acid_halide'],
                    'polyamide_C':['secondary_amide', 'NC=O-heterocycle'],
                    'polyester':['ester', 'aliphatic_alcohol', 'carboxilic_acid'],
                    'polyester_C':['eter', 'ester', 'carboxilic_acid', 'O-heterocycle'],
                    'polyurethane':['aliphatic_alcohol', 'cyanate'],
                    'polyether':['eter', 'O-heterocycle']
            }
            classes_final = []
            classes = []
            functions = []
            for i in id_number:
                for atom_mapping, o_f in zip(atom_map, org_functions):
                    if i in atom_mapping:
                        if o_f == 'alkene':
                            o_f_ = 'alkene'
                            functions.append(o_f_)
                        elif o_f == 'alkyne':
                            o_f_ = 'alkyne'
                            functions.append(o_f_)
                        elif o_f == 'primary_amine':
                            o_f_ = 'primary_amine'
                            functions.append(o_f_)
                        elif o_f == 'acid_halide':
                            o_f_ = 'acid_halide'
                            functions.append(o_f_)
                        elif o_f == 'aliphatic_alkyl_halide':
                            o_f_ = 'aliphatic_alkyl_halide'
                            functions.append(o_f_)
                        elif o_f == 'aliphatic_alcohol':
                            o_f_ = 'aliphatic_alcohol'
                            functions.append(o_f_)
                        elif o_f == 'carboxilic_acid':
                            o_f_ = 'carboxilic_acid'
                            functions.append(o_f_)
                        elif o_f == 'secondary_amide':
                            o_f_ = 'secondary_amide'
                            functions.append(o_f_)
                        elif o_f == 'NC=O-heterocycle':
                            o_f_ = 'NC=O-heterocycle'
                            functions.append(o_f_)                   
                        elif o_f == 'cyanate':
                            o_f_ = 'cyanate'
                            functions.append(o_f_)  
                        elif o_f == 'eter':
                            o_f_ = 'eter'
                            functions.append(o_f_)  
                        elif o_f == 'ester':
                            o_f_ = 'ester'
                            functions.append(o_f_)  
                        elif o_f == 'O-heterocycle':
                            o_f_ = 'O-heterocycle'
                            functions.append(o_f_) 
                
                        for key, value in classes_dict.items():
                            if set(value).issubset(functions):
                                if key == 'polyamide_A':
                                    classes.append('polyamide')
                                    break
                                if key == 'polyamide_C':
                                    classes.append('polyamide')
                                    break
                                if key == 'vinyl_A':
                                    classes.append('vinyl')
                                    break
                                if key == 'polyester_C':
                                    classes.append('polyester')
                                    break
                                else:
                                    classes.append(key)
                                    break
                            else:
                                continue
                if len(classes) == 1:
                    classes_final.append(classes)
                elif len(classes) >= 2:
                    classes_final.append(classes[0])

                if len(classes_final) == 1:
                    break

        print(f'classes_final {classes_final}')
        return classes_final

    def _open_ring(self, classes, smiles):
        '''
        Non-public method. Open ring of cyclic monomers.
        
        Arguments:
            smiles(mol object): smiles obtained as mol object from monomer.
            classes(str): classes obtained from df.
        '''
        smiles_open = []
        if classes == 'polyamide':
            if bool(re.search('(^C)1|(?<=O=C)1|(?<=N)1|(?<=\(\[:\*h\]\))1|(?<=\(\[:\*t\]\))1', str(smiles))) == True:
                smiles_open_ring = re.sub('(^C)1|(?<=O=C)1|(?<=N)1|(?<=\(\[:\*h\]\))1|(?<=\(\[:\*t\]\))1', '', str(smiles))
                smiles_open.append(smiles_open_ring)
            else:
                smiles_open.append('error03: could not open ring')
        elif classes == 'polyether':
            if bool(re.search('^C1|O1$|C1$|^O1', str(smiles))) == True:
                smiles_open_ring = re.sub('(?<=^C)1|(?<=O\(\[:\*h\]\))1$|(?<=C\(\[:\*h\]\))1$|(?<=O\(\[:\*t\]\))1$|(?<=C\(\[:\*t\]\))1$|(?<=^O)1', '', str(smiles))
                smiles_open.append(smiles_open_ring)
            elif bool(re.search('O1|C1|C1$|O1|O1$|C1', str(smiles))) == True:
                smiles_open_ring1 = re.sub('(?<=C)1(?=\w)|(?<=O)1(?=\w)', '(', str(smiles))
                smiles_open_ring2 = re.sub('(?<=O\(\[:\*h\]\))1$|(?<=C\(\[:\*h\]\))1$|(?<=O\(\[:\*t\]\))1$|(?<=C\(\[:\*t\]\))1$', ')', str(smiles_open_ring1))
                smiles_open.append(smiles_open_ring2)
            elif bool(re.search('\(\[:\*h\]\)1|\(\[:\*t\]\)1', str(smiles))) == True:
                smiles_open_ring = re.sub('(?<=\(\[:\*h\]\))1|(?<=\(\[:\*t\]\))1', '', str(smiles))
                smiles_open.append(smiles_open_ring)
            else:
                smiles_open.append('error03: could not open ring')
        elif classes == 'polyester':

            if bool(re.search('(?<=CC)\d|(?<=O)\d|(?<=\(\[:\*h\]\))\d|(?<=\(\[:\*t\]\))\d', str(smiles))) == True:
                smiles_open_ring = re.sub('(?<=CC)\d|(?<=O)\d|(?<=\(\[:\*h\]\))\d|(?<=\(\[:\*t\]\))\d', '', str(smiles))
                smiles_open.append(smiles_open_ring)
            elif bool(re.search('2|1', str(smiles))) == True:
                ester = re.search('(?=C1).*(?<=\(\[:\*h\]\))', str(smiles))
                try:
                    ester_par = '(' + ester[0] + ')'
                    ester_par_san = re.sub('1', '', str(ester_par))
                    del_ester = re.sub('(?=C1).*(?<=\(\[:\*h\]\))', '', str(smiles))
                    new_smiles = re.sub('1', ester_par_san, str(del_ester))
                    new_smiles_sanitized = re.sub('2', '1', str(new_smiles))
                    smiles_open.append(new_smiles_sanitized)
                except TypeError:
                    smiles_open.append('error04: Nonetype could not be analyzed')

            else:
                smiles_open.append('error03: could not open ring')
        return smiles_open

    def _sanitize_polymer(self, classes, smiles):
        '''
        Non-public method. Change the smiles format after opening ring.
        
        Arguments:
            smiles(mol object): smiles obtained as mol object from monomer.
            classes(str): classes obtained from df.
        '''
        smiles_sanitized = []
        if classes == 'polyamide':
            if bool(re.search('(?<=^C)1', str(smiles))) == True:
                end_carbon = re.search('(?<=\(\[:\*h\]\))(.*)', str(smiles))
                del_carbon = re.sub('(?<=\(\[:\*h\]\))(.*)|(?<=^C)1', '', str(smiles))
                new_smiles_one = end_carbon[0] + del_carbon
                new_smiles = re.sub('=O', '(=O)', str(new_smiles_one))
                smiles_sanitized.append(new_smiles)
            else:
                smiles_sanitized.append('error02: could not sanitize')

        elif classes == 'polyether':
            if bool(re.search('(?<=^C)1|(?<=^\[C\])1', str(smiles))) == True:
                end_carbon = re.search('(?<=\(\[:\*h\]\))(.*)', str(smiles))
                del_carbon = re.sub('(?<=\(\[:\*h\]\))(.*)|(?<=^C)1|(?<=^\[C\])1', '', str(smiles))
                new_smiles = end_carbon[0] + del_carbon
                smiles_sanitized.append(new_smiles)
            else:
                smiles_sanitized.append('error02: could not sanitize')

        sanitized_reduce = reduce(lambda x, y: x+y, smiles_sanitized)

        return sanitized_reduce

    def _get_atom_mapping(self, smiles, atom_mapping, org_functions, classes):
        '''
        Non-public method. Gets the classification of the polymerization mechanism for more than one functional group.
        
        Arguments:
            smiles(mol object): smiles obtained as mol object from monomer.
            org_functions(list): organic functions obtained on the get_class method.
            atom_map(list): atom mappings obtained on the get_class method.
            classes(str): classes obtained from df.
        '''
        head_tail = []
        if classes == 'vinyl':
            for atom_map, o_f in zip(atom_mapping, org_functions):
                if o_f == 'alkene':
                    try:
                        head = re.sub(f'\:{atom_map[0]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\*\d{1,2}]|\))\])|H\d|H', '', str(tail)) # get only ] near head and tail tag
                        clean_vinyl = re.sub('(?:^|(?<=\(\[:\*h\]\))=)|(?:^|(?<=O\))=)|(?:^|(?<=C\))=)', '', str(cleanSmiles)) 
                        head_tail.append(clean_vinyl)
                    except IndexError:
                        head_tail.append('error01: Could not recognize mapping')
                        break
                if o_f == 'alkyne':
                    try:
                        head = re.sub(f'\:{atom_map[0]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\*\d{1,2}]|\))\])|H\d|H', '', str(tail)) # get only ] near head and tail tag
                        clean_vinyl = re.sub('#', '', str(cleanSmiles)) 
                        head_tail.append(clean_vinyl)

                    except IndexError:
                        head_tail.append('error01: Could not recognize mapping')
                        break
                    
        if classes == 'polyamide':
            for atom_map, o_f in zip(atom_mapping, org_functions):
                try:
                    if len(org_functions) == 1:
                        if o_f == 'primary_amine':
                            head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                            tail = re.sub(f'\:{atom_map[3]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            head_tail.append(cleanSmiles)
                            break
                    elif len(org_functions) >= 2:    
                        if o_f == 'primary_amine':
                            if org_functions[-1] == 'primary_amine':
                                head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                                tail = re.sub(f'\:{atom_map[3]+1}(?!\d)', '([:*t])', str(head))
                                cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                                head_tail.append(cleanSmiles)
                                break
                            else:                                                
                                head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                        elif org_functions[0]== 'aliphatic_alcohol':
                            if (o_f == 'aliphatic_alcohol'):# & (head == ""):
                                head = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*h])', str(smiles))
                        if o_f == 'acid_halide':
                            tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            clean_carbonyl = re.sub('(?<=\(=O\))Cl|\(Cl\)|(?<=\(=O\))Cl|Cl(?=C)', '', str(cleanSmiles))
                            head_tail.append(clean_carbonyl)
                        elif o_f == 'carboxilic_acid':                            
                            tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            clean_carbonyl = re.sub('(^O)|(?:^|(?<=\(=O\))O)', '', str(cleanSmiles))
                            head_tail.append(clean_carbonyl)
                        if o_f == 'aliphatic_alkyl_halide':
                            head = re.sub(f'\:{atom_map[0]+1}(?!\d)', '([:*h])', str(smiles))
                            tail = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            clean_carbonyl = re.sub('(?<=\(=O\))Cl|\(Cl\)|(?<=\(=O\))Cl|Cl(?=C)', '', str(cleanSmiles))
                            head_tail.append(clean_carbonyl)
                            break
                except (IndexError, UnboundLocalError) as error:
                        print(f'ERROR: {error}')
                        head_tail.append('error01: Could not recognize mapping')
                        break

                if  o_f == 'NC=O-heterocycle': #for heterocycles
                    try:
                        head = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                        open_smiles = self._open_ring(classes, cleanSmiles)
                        open_smiles1 = reduce(lambda x, y: x+y, open_smiles)
                        head_tail.append(open_smiles1)

                    except IndexError:
                        head_tail.append('error01: Could not recognize mapping')
                        break

        if classes == 'polyester':
            for atom_map, o_f in zip(atom_mapping, org_functions):
                try:
                    if len(org_functions) == 1:
                        if o_f == 'aliphatic_alcohol':
                            head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                            tail = re.sub(f'\:{atom_map[3]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            head_tail.append(cleanSmiles)
                            break
                    elif len(org_functions) >= 2:
                        if o_f == 'aliphatic_alcohol':
                            head = re.sub(f'\:{atom_map[3]+1}(?!\d)', '([:*h])', str(smiles))
                        if o_f == 'carboxilic_acid':
                            if 'O-heterocycle' in org_functions:
                                head = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*h])', str(smiles))
                                tail = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*t])', str(head))
                                cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                                open_smiles = self._open_ring(classes, cleanSmiles)
                                open_smiles1 = reduce(lambda x, y: x+y, open_smiles)
                                head_tail.append(open_smiles1)
                            else:
                                tail = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*t])', str(head))
                                cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                                clean_carbonyl = re.sub('(^O)|(?:^|(?<=\(=O\))O)', '', str(cleanSmiles))
                                head_tail.append(clean_carbonyl)
                        if o_f == 'aromatic_alcohol':
                            head = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*h])', str(smiles))
                            tail = re.sub(f'\:{atom_map[8]+1}(?!\d)', '([:*t])', str(head))
                            cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                            head_tail.append(cleanSmiles)
                except IndexError:
                    head_tail.append('error01: Could not recognize mapping')
                    break

        if classes == 'polyurethane':
            for atom_map, o_f in zip(atom_mapping, org_functions):
                try:
                    if o_f == 'aliphatic_alcohol':
                        head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[3]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                        head_tail.append(cleanSmiles)
                    elif o_f == 'cyanate':
                        head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[4]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                        clean_double_bond = re.sub('=(?=N)|(?<=N)=', '', str(cleanSmiles))
                        head_tail.append(clean_double_bond)
                except (IndexError, UnboundLocalError):
                        head_tail.append('error01: Could not recognize mapping')
                        break

        if classes == 'polyether':
            for atom_map, o_f in zip(atom_mapping, org_functions):
                try:
                    print(o_f)
                    if o_f == 'O-heterocycle':        
                        head = re.sub(f'\:{atom_map[1]+1}(?!\d)', '([:*h])', str(smiles))
                        tail = re.sub(f'\:{atom_map[2]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                        open_smiles = self._open_ring(classes, cleanSmiles) #check if all polyethers are going through this step
                        if bool(re.search('(?<=^C)1', str(open_smiles))) == True:
                            sanitize = self._sanitize_polymer(classes, open_smiles)
                            head_tail.append(sanitize)
                        else:
                            head_tail.append(open_smiles)
                    if o_f == 'carbonyl':     
                        print(smiles)   
                        head = re.sub(f'CH\d\:{atom_map[0]+1}(?!\d)', r'([:*h])\g<0>', str(smiles))
                        tail = re.sub(f'\:{atom_map[0]+1}(?!\d)', '([:*t])', str(head))
                        cleanSmiles = re.sub('\:\d{1,2}|\[(?=[^\:\*])|(?:^|(?<=[\]|\d{1,2}]|\))\])|H\d|H', '', str(tail))
                        head_tail.append(cleanSmiles)
                except IndexError:
                        head_tail.append('error01: Could not recognize mapping')
                        break

        return head_tail