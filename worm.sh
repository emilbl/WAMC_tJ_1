#!/bin/bash

####
#### read input arguments
####

# variables
parameter_file=""

# usage printer function
print_usage() {
  printf "Usage:  -i <string>   Path to parameter file (required) \n"
  printf "        -h            Prints this menu \n"
}

# read input arguments
while getopts 'i:c:qph' flag; do
  case "${flag}" in
    i) parameter_file=${OPTARG} ;;
    h) print_usage
       exit 1 ;;
    *) print_usage
       exit 1 ;;
  esac
done


####
#### ensure a parameter file path is provided
####
if [ -z $parameter_file ] ; then
  echo "worm.sh: no input file supplied  ->  EXIT"
  exit 0
fi


####
#### create output directory and copy the arguments file
####
output=$(python3 create-output-dir.py $parameter_file $prev_job_folder 3>&1 1>&2 2>&3 3>&-)
echo -----[output from python]-----
printf '%s\n' "$output"
echo ------------------------------


####
#### Catch potential errors
####
if [[ $output == *"error"* ]] || [[ $output == *"Error"* ]]; then
  echo "worm.sh: the python script exited with an error  ->  EXIT"
  exit 0
fi


####
#### everything seems to be working
####
numJobs=$(echo "$output" | wc -l)
jobs=$output


####
#### start the job
####
for job in ${jobs[@]}
do
  splitted=(${job//;/ })
  job_id=${splitted[1]}
  job_path=${splitted[0]}

  bin/worm $job_path $job_id
done
