import subprocess
import sys
import pathlib
import os
import logging
import pandas as pd
import stat

EXECUTABLE_PATH = ".\\x64\\Release\\TSP_Project.exe"

METHODS = {
    "hardfix_tuning" : ["-solver hf -hf2opt 1 -hfpstart 0.4 -hfpscal 0.03", "-solver hf -hf2opt 1 -hfpstart 0.4 -hfpscal 0.00", "-solver hf -hf2opt 1 -hfpstart 0.4 -hfpscal 0.02","-solver hf -hf2opt 1 -hfpstart 0.4 -hfpscal 0.05", "-solver hf -hf2opt 1 -hfpstart 0.6 -hfpscal 0.05", "-solver hf -hf2opt 1 -hfpstart 0.25 -hfpscal 0.02"],
    "hardfix_tuning_2" : ["-solver hf -hfpstart 0.5 -hfpscal 0.00", "-solver hf -hfpstart 0.5 -hfpscal 0.03", "-solver hf -hfpstart 0.5 -hfpscal 0.06","-solver hf -hfpstart 0.5 -hfpscal 0.10", "-solver hf -hfpstart 0.3 -hfpscal 0.00", "-solver hf -hfpstart 0.3 -hfpscal 0.03", "-solver hf -hfpstart 0.3 -hfpscal 0.06", "-solver hf -hfpstart 0.3 -hfpscal 0.10"],
    "metaheuristics" : ["-solver vns -vns_min_kicks 3 -vns_max_kicks 3", "-solver tabu -varten 1 -tabuavg 3.0"],
    "branch_and_cut" : ["-solver bc", "-solver bcf", "-solver bcm", "-solver bcfm"],
    "tabu" : ["-solver tabu -varten 0 -tabuavg 1.0" , " -solver tabu -varten 0 -tabuavg 3.0 "," -solver tabu -varten 0 -tabuavg 4.0 " , "-solver tabu -varten 0 -tabuavg 2.0" , " -solver tabu -varten 1 -tabuavg 2.0", "-solver tabu -varten 1 -tabuavg 3.0"],
    "heuristics" : [ "-solver nn -greperc 0.0", "-solver nn -greperc 1.0", "-solver 2opt -greperc 0.0333" , "-solver 2opt -greperc 1.0"],
    "tabuavg" : [ "-solver tabu -tabuavg 1.5" ,"-solver tabu -tabuavg 2.0" , "-solver tabu -tabuavg 2.5", "-solver tabu -tabuavg 3" , "-solver tabu -tabuavg 4" ],
    "tabuamp" : [ "-solver tabu -tabuamp 0.05" ,"-solver tabu -tabuamp 0.1" , "-solver tabu -tabuamp 0.2", "-solver tabu -tabuamp 0.3"],
    "tabufreq" : [ "-solver tabu -tabufreq 0.1" , "-solver tabu -tabufreq 1.0"],
    #5 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14
    "exact" : ["-solver benders", "-solver bc" , "-solver bcm" , " -solver bcp" , " -solver bcmp" , " -solver bcf" , " -solver bcfm" , " -solver bcfp" , " -solver bcfmp"],
    "math" : ["-solver hf -hfpstart 0.5 -hfpscal 0.10" , "-solver sf -sfkstart 20 -sfkscal 15"],
    "vns_tune" : ["-solver vns -vns_min_kicks 1 -vns_max_kicks 1", "-solver vns -vns_min_kicks 3 -vns_max_kicks 3","-solver vns -vns_min_kicks 5 -vns_max_kicks 5", "-solver vns -vns_min_kicks 1 -vns_max_kicks 3", "-solver vns -vns_min_kicks 1 -vns_max_kicks 5"],
    "softfixing_tuning" : ["-solver sf -sfkstart 20 -sfkscal 0" , "-solver sf -sfkstart 20 -sfkscal 5" ,"-solver sf -sfkstart 20 -sfkscal 10" ,"-solver sf -sfkstart 20 -sfkscal 15" ,"-solver sf -sfkstart 20 -sfkscal 20"]
}

TEST_BED_SIZE = 15

RANDOM_SEED = 0

NNODES = 700

TIMELIMIT = 60

PERFPLOT_TYPE = ["cost", "time"] #do not touch

VERBOSITY = 0

def generate_csv_col_name(type)-> str:
    # ciclo tipo parse cmd line 
    col_name =""
    first_item = 1
    tokens = type.split()
    for token in tokens:
        if (token.startswith("-")):
            continue
        else:
            if first_item:
                col_name = col_name + f"{token}"
                first_item =0
            else:
                col_name = col_name + f"_{token}"
    
    return col_name
    
if __name__ == "__main__":
    
    if not (os.path.exists(".\\perfplot\\")):
        os.makedirs(".\\perfplot\\")
    
    if not (os.path.exists(".\\runs\\")):
        os.makedirs(".\\runs\\")

    logger = logging.getLogger('Travelling Salesman Problem')
    logger.setLevel(logging.INFO)
    # Parse input to pyhton file 
    if len(sys.argv) < 2:
        print("Usage: python script.py <method>")
        sys.exit()
    
    # parse for method in command line
    type = sys.argv[1]
    if type not in METHODS.keys():
        print("method not recognized")
        sys.exit()
    
    csv_column_names =[]
    
    csv_instance_id = f"{NNODES}_{RANDOM_SEED}_{TEST_BED_SIZE}"

    
    for tsp_run in METHODS[type]:
        # generate single column file 

        csv_column_names.append(generate_csv_col_name(tsp_run))

        # run tsp instance 
        try:
            for pp_type in PERFPLOT_TYPE:                
                file = open(f".\\runs\\{csv_column_names[-1]}_{csv_instance_id}_{pp_type}.csv", 'w')
                file.close()
                
            print(csv_column_names[-1])
            str_exec = f"{EXECUTABLE_PATH} -nnodes {NNODES} -seed {RANDOM_SEED} -tl {TIMELIMIT} -verbose {VERBOSITY} {tsp_run} -csv_column_name {csv_column_names[-1]} -test_bed_size {TEST_BED_SIZE}"
            print("Running " + tsp_run)
            output = subprocess.run((str_exec), capture_output=True, text=True).stdout

            
            
            #os.system(str_exec)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing {tsp_run} : {e}")
            continue
        except IndexError as e:
            logger.error(f"Index Error {e}")
            continue

    # Now I have the single column files: for time and cost we initializa new df 
    to_be_perfplotted = []
    for pp_type in PERFPLOT_TYPE:

        pp_df = pd.DataFrame(columns = ([f'{len(METHODS[type])}'] + csv_column_names))
        print("columns : "+ pp_df.columns)
        pp_csv_name = f'{pp_type}_{type}_{csv_instance_id}.csv'
        
        # Create the df for perfplot 
        for column in csv_column_names:
            df = pd.read_csv((f'.\\runs\\{column}_{csv_instance_id}_{pp_type}.csv'), sep=', ')
            print(df)
            pp_df[column] = df[column]
        
        #first column init 
        pp_df[f'{len(METHODS[type])}'] = df['1']
    
        #TODO: pp_Df to csv
        path = f".\\perfplot\\{pp_csv_name}"
        pp_df.to_csv(pathlib.Path(path), index = False)
        to_be_perfplotted.append(path)
    
    #generate perf plot 
    for name in to_be_perfplotted:
        pp_name =  name.replace(".csv", ".pdf")
        if ("time" in pp_name):
            str_exec = f"python .\\perfprof.py -D , -M 8.0 -X \"Time Ratio\" {name} {pp_name}"
        else:
            str_exec = f"python .\\perfprof.py -D , -M 1.1 -X \"Cost Ratio\" {name} {pp_name}"
        os.system(str_exec)