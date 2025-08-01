import sys
import os
from datetime import datetime
from subprocess import call


# HERE GIVE THE FULL PATH TO GAMESS FOLDER:
path_to_gamess = "/dccstor/rgiro1/gamess" #DEFAULT

# HERE GIVE THE VERSION OF GAMESS YOU HAVE:
version = "00" 

path_input_file = sys.argv[1]
print('path_input_file',path_input_file)
	
try:
	number_of_processors = sys.argv[2] 
except:
	number_of_processors = 1

input_directory = os.getcwd()

output_file = input_directory + "/" + path_input_file.rsplit("/", 1)[0]
print('output_file ',output_file)
name_file = path_input_file.split("/")[2]
name = name_file.split(".inp")[0]
output_name = name+".log"

#call(["cp", path_input_file, path_to_gamess], shell=True)

#os.chdir(path_to_gamess)
	
# RUNNING GAMESS JOB
start_time = datetime.now()
#call(["rungms", name_file, version, number_of_processors, ">", output_name], shell=True)
os.system(path_to_gamess + '/rungms ' + output_file +'/'+ name_file +' '+ version + ' 1 > ' + output_name)

end_time = datetime.now()
print(f'GAMESS simulation for {name_file} took {end_time - start_time}')

#call(["copy", output_name, output_file], shell=True)
os.system('cp '+ output_name + ' '+ output_file)

# CHECKING THAT THERE ARE NO OLD RESIDUAL FILES
#try:
#	call(["del", path_to_gamess+"\\"+name+".inp"], shell=True)
#	call(["del", path_to_gamess+"\\"+name+".log"], shell=True)
#except:
#	print(">> No residual files found.")
os.system('rm '+output_name)

print(f"GAMESS simulation for {name_file} is completed.")
