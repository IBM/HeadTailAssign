import pandas as pd
import csv
import re


class Helper:
    '''Class related to helper functions.

    Methods:
        _get_atoms_list(self, path): Non-public method. Gets the atom type on the file .txt. Returns a list.
        _get_n_orbitals(self, path): Non-public method. Gets the number of orbitals calculated by GAMESS. Returns a list.
        _get_id(self, path): Non-public method. Gets the atom id. Returns a list.
        _get_a_population(self, path, n_orbitals): Non-public method. Gets the atomic population calculated by GAMESS. Returns a list.
        _change_atom_ID(self, atomId): Changes the identification number of the atoms by minus one. Returns a list
        _get_dataframe(self, valence_orbital_list, filename): Gets a dataframe from lists with informations about the atoms. Returns a dataframe.
        _merge_head_tail(self, df): Merges the input dataframe with the output dataframe. Returns a dataframe.
    '''

    def _get_atoms_list(self, path):
        '''Non-public method. Gets the atom type on the file .txt. Returns a list.
        
        Arguments:
            path: path file with the information regarding atom types.
            '''
        with open(path, encoding='utf8') as f:
            should_append = False
            atoms = []
            for x, line in enumerate(f):
                if line.find('          INTERNUCLEAR DISTANCES (ANGS.)')!=-1:
                    should_append=False
                if should_append:
                    if line[1:3]:
                        atoms.append(line[1:3])
                if line.find('     ATOMS')!=-1:
                    should_append=True
            atoms_clean = [a.strip(' ') for a in atoms]
        return atoms_clean
    
    def _get_n_orbitals(self, path):
        '''Non-public method. Gets the number of orbitals calculated by GAMESS. Returns a list.
        
        Arguments:
            path: path file with the information regarding orbitals.
            '''
        with open(path, encoding='utf8') as f:
            n_orbitals = []
            for x, line in enumerate(f):
                if re.search('ORBITALS ARE OCCUPIED', line):
                    n_orbitals.append(line[4:7])
            n_orbitals_clean = [o.strip(' ') for o in n_orbitals]
            
        return n_orbitals_clean
    
    def _get_id(self, path):
        '''Non-public method. Gets the atom id. Returns a list.
        
        Arguments:
            path: path file with the information of atomic id.
            '''
        with open(path, encoding='utf8') as f:
            n_id = []
            should_append = False
            for x, line in enumerate(f):
                if should_append:
                    if line:
                        n_id.append(line[0:5])
                    else:
                        break
                if line.find('---------------------     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL')!=-1:
                    should_append=True
            n_id_clean = [i.strip(' ') for i in n_id]
            n_id_clean_a = []
            for n in n_id_clean:
                if n == '\n':
                    continue
                if n == '':
                    continue
                else:
                    if n not in n_id_clean_a:
                        n_id_clean_a.append(n)

        return n_id_clean_a
    
    def _get_a_population(self, path, n_orbitals, n_id):
        '''Non-public method. Gets the atomic population calculated by GAMESS. Returns a list.
        
        Arguments:
            path: path file with the information atomic population. UPDATE HERE
            '''
        with open(path, encoding='utf8') as f:
            o_p = []
            should_append = False
            count = 0
            div_rest = (int(n_orbitals[0]) % 5)
            n_id_compare = []
            for x, line in enumerate(f):
                
                if should_append:
                    for n in n_id:
                        o_id = line[0:6]
                        c_line = o_id.replace(" ", "")
                        if n == c_line:
                            n_id_compare.append(c_line)
                        elif n_id_compare == n_id:
                            count += 1
                            n_id_compare = []
                        if count == div_rest -1:
                            if n == c_line:
                                if div_rest == 1:
                                    o_p.append(line[18:27].replace("\n", ""))
                                if div_rest == 2:
                                    o_p.append(line[29:38].replace("\n", ""))
                                if div_rest == 3:
                                    o_p.append(line[40:49].replace("\n", ""))
                                if div_rest == 4:
                                    o_p.append(line[51:60].replace("\n", ""))
                                if div_rest == 5:
                                    o_p.append(line[62:71].replace("\n", ""))
                if line.find('---------------------     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL')!=-1:
                    should_append=True
            
        return o_p

    def _change_atom_ID(self, atomId):
        '''Non-public method. Changes the identification number of the atoms by minus one. Returns a list.
        
        Arguments:
            atomId (list): Identification number of the atoms.
            ''' 
        atomIdList = []
        for i in atomId:
            newID = int(i)
            atomIdList.append(newID)
        
        return atomIdList

    def _get_dataframe(self, atoms, n_id, a_population, filename):
        '''Non-public method. Gets a dataframe from lists with informations about the atoms. Returns a dataframe.
        
        Arguments:
            valence_orbital_list (list): List with only valence orbitals.
            filename (str): Name of the files without the extension.
            ''' 

        # atomIdList = self._change_atom_ID(atomId)
       
        df = pd.DataFrame(
            {'atom': atoms,
            'id': n_id,
            'Mulliken population': a_population}
        )
        df.sort_values(by=['Mulliken population'], ascending=False, inplace=True)
        df.to_csv(f"{filename}.csv", sep=",", index=False, quoting=csv.QUOTE_ALL)
        
        return df

    def _merge_head_tail(self, df):
        '''Non-public method. Merges the input dataframe with the output dataframe. Returns a dataframe.
        
        Arguments:
            df (dataframe): Input dataframe. 
            ''' 
        df1 = pd.read_csv('input.csv', sep=',', names=["name","reaction","monomers","classes","mechanism"])
        df2 = df1.merge(df, on='name')
        df2.to_csv(f"head_tail_comparison.csv", sep=",", index=False, quoting=csv.QUOTE_ALL)

        return df2