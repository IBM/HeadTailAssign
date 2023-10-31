from collections import defaultdict
import glob, os
import numpy as np
import re
import pandas as pd
import csv
import shutil

from HeadTailAssign.extractor_helper import Helper


class Extractor:
    '''Class related to data extraction.

    Methods:
        op_extractor(self): Extracts only the data regarding Mulliken Atomic Overlap Populations from .out file. Saves .txt file. Returns None.
        rx_extractor(self): Extracts the Rx value organized by higher to lower value. Returns a dataframe.
        get_head_tail(self): Gets head and tail data from a .csv file and creates a dataframe. Returns a dataframe.
        analyze_head_tail(self): Analyzes the head and tail data to compare if the information is equal (add True) or different (add False). Returns a dataframe.
        '''
    def op_extractor(self, name_dir: str): 
        '''Extracts only the data regarding Mulliken Atomic Overlap Populations from .out file. Saves .txt file. Returns None.
         
        Arguments:
            name_dir(str) = Name of the output directory'''
        try:
            for root, dirs, files in os.walk(name_dir):
                for file in files:
                    if file.endswith(".log"):
                        fileName = re.sub(".log", "", file)
                        path = root + "\\" + file 
                        path_output = root + "\\" + fileName + ".txt"

                        #https://stackoverflow.com/questions/61172255/how-to-slice-data-from-a-text-file-given-the-desired-range-of-lines
                        with open(path, encoding='utf8') as f:
                            with open(path_output, "w") as out:
                                print(f'Analyzing {file}...')
                                out.write(f'{file}\n')
                                should_write=False
                                for x, line in enumerate(f):
                                    if line.find('               ----- POPULATIONS IN EACH AO -----')!=-1:
                                        should_write=False
                                    if should_write:
                                        out.write(line)
                                    if line.find('     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL')!=-1:
                                        out.write('     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL\n')
                                        should_write=True
                                    
                                    if line.find(' SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).')!=-1:
                                        out.write('\n     OCCUPIED ORBITALS\n')
                                        for i in range(1):
                                            out.write(next(f))
                                            out.write(f'\n---------------------')
                                    
                                    if line.find('          INTERNUCLEAR DISTANCES (ANGS.)')!=-1:
                                        should_write=False
                                    if line.find('           CHARGE         X                   Y                   Z')!=-1:
                                        out.write('\n     ATOMS\n')
                                        should_write=True

                                        
                            if os.stat(path_output).st_size >= 100:
                                print(f'file {fileName}.txt was generated sucessfully.')
                            else:
                                f.close()
                                print(f'WARNING: {fileName} could not be generated.')
                                os.remove(path_output)
                                                
        except ValueError:
            raise ValueError("There is no files with .log extension.")

        return
    # def op_extractor(self, name_dir: str): 
    #     '''Extracts only the data regarding Mulliken Atomic Overlap Populations from .out file. Saves .txt file. Returns None.
         
    #     Arguments:
    #         name_dir(str) = Name of the output directory'''
    #     try:
    #         for root, dirs, files in os.walk(name_dir):
    #             for file in files:
    #                 if file.endswith(".log"):
    #                     fileName = re.sub(".log", "", file)
    #                     path = root + "\\" + file 
    #                     path_output = root + "\\" + fileName + ".txt"

    #                     #https://stackoverflow.com/questions/61172255/how-to-slice-data-from-a-text-file-given-the-desired-range-of-lines
    #                     with open(path, encoding='utf8') as f:
    #                         with open(path_output, "w") as out:
    #                             print(f'Analyzing {file}...')
    #                             out.write(f'{file}\n')
    #                             out.write('          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----\n')
    #                             should_write=False
    #                             for x, line in enumerate(f):
    #                                 if line.find('          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----')!=-1:
    #                                     should_write=False
    #                                 if should_write:
    #                                     out.write(line)
    #                                 if line.find('               ----- POPULATIONS IN EACH AO -----')!=-1:
    #                                     should_write=True

    #                         if os.stat(path_output).st_size >= 100:
    #                             print(f'file {fileName}.txt was generated sucessfully.')
    #                         else:
    #                             f.close()
    #                             print(f'WARNING: {fileName} could not be generated.')
    #                             os.remove(path_output)
                                                
    #     except ValueError:
    #         raise ValueError("There is no files with .log extension.")

    #     return

    def rx_extractor(self, name_dir):
        '''Extracts the Rx value organized by higher to lower value. Returns a dataframe.
        
        Arguments:
            name_dir(str) = Name of the output directory.'''
        global results
        try:
            for root, dirs, files in os.walk(name_dir):
                for file in files:
                    if file.endswith(".txt"):
                        fileName = re.sub(".txt", "", file)
                        path = root + "\\" + file 
                        path_output = root + "\\" + fileName

                        if os.stat(path).st_size >= 100:
                            print(f'Analyzing {file}...')
                            #https://stackoverflow.com/questions/33538660/how-do-i-compare-2-lines-in-a-string-in-python
                            with open(path, encoding='utf8') as f:
                                valence_list = []
                                for x, line in enumerate(f):
                                    if re.search(r'\W \d', line):
                                        values = line[17:42]

                                        if bool(re.search(r'[^XYZ]Z', values[5:7])) == True:
                                            z_orbital = values
                                            valence_list.append(z_orbital)

                            helper = Helper()
                            if len(valence_list) <= 1:
                                ("WARNING: There is no atoms with z orbital in this monomer. Files will be deleted.")
                                shutil.rmtree(root)
                                break
                            else:
                                valence_orbital_list = helper._get_valence_orbital(valence_list)
                                valence_orbital_list_updated = helper._get_valence_orbital_last(valence_list, valence_orbital_list)
                                results = helper._get_dataframe(valence_orbital_list_updated, path_output)

        except ValueError:
            raise ValueError("There is no files with .txt extension.")

        return results

    def get_head_tail(self, name_dir: str):
        '''Gets head and tail data from a .csv file and creates a dataframe. Returns a dataframe.

        Arguments:
            name_dir(str) = Name of the output directory.'''
        os.chdir(name_dir)
        d = defaultdict(list)
        for file in glob.glob("*.csv"):
            if file == 'input.csv':
                pass
            else:
                fileName = re.sub(".csv", "", file)
                df = pd.read_csv(file)
                d['name'].append(fileName)
                d['head output'].append(df.iloc[1,1])
                d['tail output'].append(df.iloc[0,1])
        df1 = pd.DataFrame.from_dict(d)

        helper = Helper()
        helper._merge_head_tail(df1)

        os.chdir("../")
        return df1

    def analyze_head_tail(self):
        '''Analyzes the head and tail data to compare if the information is equal (add True) or different (add False). Returns a dataframe.'''
        df = pd.read_csv('head_tail_comparison.csv', sep=',')
        #https://stackoverflow.com/questions/27474921/compare-two-columns-using-pandas
        df['comparison'] = np.where((df['head input'] == df['head output']) & (df['tail input'] == df['tail output'])
                     , 'True', 'False')
        df.to_csv(f"smi2gamess-output-comparison.csv", sep=",", index=False, quoting=csv.QUOTE_ALL)

        return df