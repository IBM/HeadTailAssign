import pandas as pd
import re
from IPython.display import display
from rdkit import Chem
import networkx as nx


class Performance:
    '''
    Class related to data performance calculation

    Methods:

        _compare(self, df_merged) -> (list): Private method. Compare SMILES strings. Returns list.
        _canonicalize(self, df_merged): Private method. Performs the canonicalization of SMILES strings. Returns Dataframe.
        validate_results(self, df1_smiles, df2_smiles): Performs the validation of the results by comparing predicted data with actual data. Returns Dataframe.
        calculate_performance_name(self): Calculates the performance of the m2p algorithm. It uses the formula data_calculated/total_data*100. Returns string.
        get_not_calculated(self): Get only the polymers which m2p could not polymerize. Returns dataframe.
        
        '''
    def __init__(self, df1, df2):
        '''Initialize the instance of a class.

        Arguments:
            df1(dataframe): Dataframe which contains m2p output.
            df2(dataframe): Dataframe which contains m2p input.
        '''
        self.df1 = df1
        self.df2 = df2

    def _topology_from_rdkit(self, rdkit_molecule):

        topology = nx.Graph()
        for atom in rdkit_molecule.GetAtoms():
            # Add the atoms as nodes
            topology.add_node(atom.GetIdx())

            # Add the bonds as edges
            for bonded in atom.GetNeighbors():
                topology.add_edge(atom.GetIdx(), bonded.GetIdx())

        return topology


    def _is_isomorphic(self, topology1, topology2):
        return nx.is_isomorphic(topology1, topology2)

    def _compare(self, df_merged) -> (list):
        '''Private method. 
        Compare SMILES strings. Returns list.
        
        Arguments:
            df_merged(dataframe): Dataframe which contains merged data from df1 and df2.'''

        print(df_merged)
        df_rows = len(df_merged)

        n = 0

        number = []
        name = []
        tf = []
        smiles_pred = []
        smiles_actual = []

        while n < df_rows:
            try:
                a = self._is_isomorphic(self._topology_from_rdkit(Chem.MolFromSmiles(df_merged['canonicalized_x'][n][0])), self._topology_from_rdkit(Chem.MolFromSmiles(df_merged['canonicalized_y'][n][0])))
                # print(df_merged['canonicalized_x'][n][0], df_merged['canonicalized_y'][n][0])
                #a = Chem.MolFromSmiles(df_merged['canonicalized_x'][n][0]).HasSubstructMatch(Chem.MolFromSmiles(df_merged['canonicalized_y'][n][0]))
                #a = set(df_merged['canonicalized_x'][n]).issubset(df_merged['canonicalized_y'][n])
                name.append(df_merged['name'][n])
                smiles_pred.append(df_merged['canonicalized_y'][n])
                smiles_actual.append(df_merged['canonicalized_x'][n])
                number.append(n)
                tf.append(a)
            except:
                name.append(df_merged['name'][n])
                smiles_pred.append(df_merged['canonicalized_y'][n])
                smiles_actual.append(df_merged['canonicalized_x'][n])
                number.append(n)
                tf.append(False)               
            n += 1

        df_res = pd.DataFrame({'name': name, 'number': number, 'smiles_pred': smiles_pred, 'smiles_actual': smiles_actual, 'boolean': tf})
        print(df_res)

        res = []
        name = []
        s_pred = []
        s_actual = []
        def f(x, y, z, a, b):
            if y == False:
                res.append(x)
                name.append(z)
                s_pred.append(a)
                s_actual.append(b)

        result = [f(x, y, z, a, b) for x, y, z, a, b in zip(df_res['number'], df_res['boolean'], df_res['name'], df_res['smiles_pred'], df_res['smiles_actual'])]
        print(f'{len(res)} results are False, those are: {name}')

        return df_res

    def _canonicalize(self, df_merged):
        '''Private method. 
        Performs the canonicalization of SMILES strings. Returns Dataframe.
        
        Arguments:
            df_merged(dataframe): Dataframe which contains merged data from df1 and df2.'''
        n = 0

        df_rows = len(df_merged)

        df1_can_res = []
        df2_can_res = []

        while n < df_rows:
            # print(df_merged['head_tail_y'][n][2:-2])
            # print(Chem.MolToSmiles(Chem.MolFromSmiles(df_merged['head_tail_y'][n][2:-2]), True))
            try:
                df1_can = Chem.MolToSmiles(Chem.MolFromSmiles(df_merged['head_tail_x'][n]),True)
                #print(df1_can)
                df1_can_res.append(df1_can)
                #print(df1_can_res)
            except:
                df1_can_res.append('no')
            try:
                df2_can = Chem.MolToSmiles(Chem.MolFromSmiles(df_merged['head_tail_y'][n][2:-2]),True)
                df2_can_res.append(df2_can)
            except:
                df2_can_res.append('no')

            n += 1

        df_merged['canonicalized_x'] = df1_can_res
        df_merged['canonicalized_y'] = df2_can_res

        return

    def validate_results(self, df1_smiles, df2_smiles):
        '''Performs the validation of the results by comparing predicted data with actual data. Returns Dataframe.
        
        Arguments:
            df1_smiles(str): Name of the SMILES column of df1.
            df2_smiles(str): Name of the SMILES column of df2.'''
        
        df1 = self.df1.sort_values([df1_smiles])
        df2 = self.df2.sort_values([df2_smiles])
        df_merged = df2.merge(df1, on=df1_smiles)

        self._canonicalize(df_merged)

        df_res = self._compare(df_merged)

        return df_res

    def calculate_performance(self) -> (str):
            '''Calculates the performance of the m2p algorithm. It uses the formula data_calculated/total_data*100. Returns string.'''

            not_classes = self.df2[(self.df2["classes"] != "not recognized") & (~self.df2.classes.str.contains('error'))]
            not_mechanism = not_classes[(not_classes["mechanism"] != "not recognized") & (~not_classes.mechanism.str.contains('error'))]
            not_head_tail = not_mechanism[(not_mechanism["head_tail"] != "not recognized") & (~not_mechanism.head_tail.str.contains('error'))]
            df1_rows = len(self.df1.index)
            df2_rows = len(not_head_tail.index)

            result = round((df2_rows / df1_rows)*100)
            print(f'Performance = {result}%')
            
            return result

    def get_not_calculated(self):
        '''Get only the polymers which m2p could not polymerize. Returns dataframe.'''

        merged = self.df1.merge(self.df2, on='head_tail', indicator=True, how='outer')
        merged = merged[merged['_merge'] != 'both']

        merged_rows = len(merged)
        print(f'm2p algorithm could not calculate {merged_rows} polymers')

        return merged