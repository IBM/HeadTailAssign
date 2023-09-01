import sys
import os
from datetime import datetime
from subprocess import call


# HERE GIVE THE FULL PATH TO GAMESS FOLDER:
path_to_gamess = "C:\\Users\\Public\\gamess-64" #DEFAULT

# HERE GIVE THE VERSION OF GAMESS YOU HAVE:
version = "2020.R2.pgiblas" 

path_input_file = sys.argv[1] 
	
try:
	number_of_processors = sys.argv[2] 
except:
	number_of_processors = 1

input_directory = os.getcwd()

output_file = input_directory + "\\" + path_input_file.rsplit("\\", 1)[0]
name_file = path_input_file.split("\\")[2]
name = name_file.split(".inp")[0]
output_name = name+".log"

call(["copy", path_input_file, path_to_gamess], shell=True)

os.chdir(path_to_gamess)
	
# RUNNING GAMESS JOB
start_time = datetime.now()
call(["rungms.bat", name_file, version, number_of_processors, ">", output_name], shell=True)

end_time = datetime.now()
print(f'GAMESS simulation for {name_file} took {end_time - start_time}')

call(["copy", output_name, output_file], shell=True)

# CHECKING THAT THERE ARE NO OLD RESIDUAL FILES
try:
	call(["del", path_to_gamess+"\\"+name+".inp"], shell=True)
	call(["del", path_to_gamess+"\\"+name+".log"], shell=True)
except:
	print(">> No residual files found.")

print(f"GAMESS simulation for {name_file} is completed.")
