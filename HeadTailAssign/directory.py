import os
import time
from art import *
import pandas as pd
import sys

class Directory:
    '''Class related to directory manipulation

    Methods
        create_output_dir(self, name: str) = Function related to directory manipulation. It will create the directory in which would be created outputfiles.
    '''

    def starting(self, disable=False):

        if disable==False:
            logo=text2art("HeadTailAssign")
            print(logo)
            print('Code developed by: Brenda Ferrari and Ronaldo Giro\n\n')
            time.sleep(1)
            print("-------------------------------------------------")
            print("Assignment is starting...")
            time.sleep(1)        

        return

    def get_data(self, args, print_df=False):

        if len(args) == 1:
            print("Some parameters are missing. Please follow the structure main.py [data]")
            sys.exit()

        data = str(args[1])
        df = pd.read_csv(data, sep=',')

        if print_df == True:
            print('\n\nThe following data is going to be analyzed:')
            print(df)
            print('')
            time.sleep(1)

        return df

    def create_output_dir(self, name: str):
        '''Function related to directory manipulation. It will create the directory in which would be created outputfiles.

        Arguments:
            name(str): name of the output directory. Default is output-directory
        '''
        cwd = os.getcwd()
        print(f'Creating output directory at {cwd}')
        if name:
            os.mkdir(name)

        return