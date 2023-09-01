import pandas as pd
import csv


class Helper:
    '''Class related to helper functions.

    Methods:
        get_valence_orbital(self, valence_list): Gets the atom valence orbital. Returns a list.
        get_valence_orbital_last(self, valence_list, valence_orbital_list): Gets the valence orbital of the last atom. Returns a list.
        change_atom_ID(self, atomId): Changes the identification number of the atoms by minus one. Returns a list
        get_dataframe(self, valence_orbital_list, filename): Gets a dataframe from lists with informations about the atoms. Returns a dataframe.
        merge_head_tail(self, df): Merges the input dataframe with the output dataframe. Returns a dataframe.
    '''

    def _get_valence_orbital(self, valence_list):
        '''Non-public method. Gets the atom valence orbital. Returns a list.
        
        Arguments:
            valence_list: List with all Z orbitals from a atom.
            '''
        valence_orbital_list = []
        for valence_orbital in range(len(valence_list)-1):

            lineOne = valence_list[valence_orbital]
            lineTwo = valence_list[valence_orbital+1]

            if (lineOne[:8] != lineTwo[:8]):
                valence_orbital_list.append(lineOne)
            elif (lineOne[:8] == lineTwo[:8]):
                valence_orbital_list.append(lineTwo)

        return valence_orbital_list

    def _get_valence_orbital_last(self, valence_list, valence_orbital_list):
        '''Non-public method. Gets the valence orbital of the last atom. Returns a list.
        
        Arguments:
            valence_list (list): List with all Z orbitals from a atom.
            valence_orbital_list (list): List with only valence orbitals.
            '''      
        for valence_orbital in range(len(valence_list)>-2):
            lastLine = valence_list[-1]
            nextLastLine = valence_list[-2]
            if (nextLastLine[:8] != lastLine[:8]):
                valence_orbital_list.append(lastLine)
                
        return valence_orbital_list

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

    def _get_dataframe(self, valence_orbital_list, filename):
        '''Non-public method. Gets a dataframe from lists with informations about the atoms. Returns a dataframe.
        
        Arguments:
            valence_orbital_list (list): List with only valence orbitals.
            filename (str): Name of the files without the extension.
            ''' 
        atomList = []
        atomId = []
        orbitalList = []
        OPValueList = []
        for valence_orbital in valence_orbital_list:
            atomList.append(valence_orbital[:2])
            atomId.append(valence_orbital[2:4])
            orbitalList.append(valence_orbital[4:8])
            OPValueList.append(valence_orbital[13:20])

        atomIdList = self._change_atom_ID(atomId)
       
        df = pd.DataFrame(
            {'atom': atomList,
            'id': atomIdList,
            'valence orbital': orbitalList,
            'Mulliken populations': OPValueList}
        )
        df.sort_values(by=['Mulliken populations'], ascending=False, inplace=True)
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