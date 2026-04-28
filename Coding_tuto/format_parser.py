import argparse
import sys
sys.path.append('../..')

#the following if condition is used such that the script can be used as a module or a main script, i.e. when it's called from another script it does not need these arguments and we can use its functions easily
if __name__ == "__main__": 
	args = parser.parse_args()
	parser = argparse.ArgumentParser(description='This script is used for...')

	parser.add_argument('--path',dest='path_var',default='./',action=,help='Provide the path of the experiment results.')

	args = parser.parse_args()
	path_experiment = args.path_var

	if args.path_var is None:
	    parser.print_help()
	    sys.exit("Error: Missing required arguments. Please provide values for --path.")

Everything that is not part of the main, i.e. of the part that would be called alone directly at the cmd line, i.e. that is outside of what
we call from another script must be in the if __main__



For example:
import argparse

# Define functions needed for both other_script.py and plot_left.py
def common_function1():
    # implementation
    pass

def common_function2():
    # implementation
    pass

# Encapsulate module-level code within an if __name__ == "__main__": block
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Description of your script')
    parser.add_argument('--path', type=str, help='Description of path argument')
    parser.add_argument('--path_db', type=str, help='Description of path_db argument')
    parser.add_argument('--variability', type=float, help='Description of variability argument')
    parser.add_argument('--bsj', type=int, help='Description of bsj argument')
    
    args = parser.parse_args()

    # Add any module-level code that depends on variables outside functions here
    path_experiment = args.path
    count_matrix = generate_count_matrix(path_experiment, False)  # Adjust this line based on your actual code

    # Continue with the rest of your script logic using args here

# Define other functions specific to other_script.py
def other_function1():
    # implementation
    pass

def other_function2():
    # implementation
    pass

# Export functions needed in plot_left.py
def check_database_circpedia_and_filter_seq_type(df_bed, value, file1, file2, flag):
    # implementation using common_function1 and common_function2 if needed
    pass

def create_df_db_circpedia(path_db):
    # implementation using common_function1 and common_function2 if needed
    pass

def intersection(list1, list2):
    # implementation using common_function1 and common_function2 if needed
    pass
