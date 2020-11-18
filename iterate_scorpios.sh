#!/bin/bash

#############################################################################################
# Bash wrapper script to run SCORPiOs in iterative mode                                     #
#                                                                                           #
# Example usage: bash iterate_scorpios.sh --snake_args="--configfile config_example.yaml"   #
#                                         [--max_iter=5] [--min_corr=1]                     #
#                                         [--starting_iter=1]                               #
#############################################################################################


#Set default values for parameters (--j and --snake_args are required)
max_iter=5
min_corr=5
iteration=1

### Command-line argument parsing ###
while [ $# -gt 0 ]; do

  case "$1" in

    #maximum number of iteration to run
    --max_iter=*)
      max_iter="${1#*=}"
      ;;

    #stop if less than min_corr trees have been corrected
    --min_corr=*)
      min_correction="${1#*=}"
      ;;

    #starting iteration, if we want to resume a SCORPiOs run for a certain iter
    --starting_iter=*)
      iteration="${1#*=}"
      ;;

    #arguments for snakemake execution
    --snake_args=*)
      snake_args="${1#*=}"
      ;;

    *)

      echo "*********************************" >&2
      echo " ArgumentError: Invalid argument:" >&2
      echo " $1                              " >&2
      echo "*********************************" >&2
      exit 1

  esac

  shift

done

#check that required arg is set
if [ -z "$snake_args" ]; then
  echo "******************************************" >&2
  echo " ArgumentError: --snake_args is required  " >&2
  echo "******************************************" >&2
  exit 1
fi

#extract config related args from snakemake args, to invoke the --config option correctly below
snake_config_args=${snake_args#*--config }

if [ "$snake_config_args" == "$snake_args" ]; then

  snake_config_args="--config "

else

  snake_config_args=${snake_config_args%% -*}
  snake_config_args="--config ${snake_config_args}"
  snake_args="${snake_args/${snake_config_args}/}"

fi

j=0

#get configfile name to then extract the jobname and the species tree from it
configfile=${snake_args#*--configfile=}

if [ "$configfile" == "$snake_args" ]; then

  configfile=${snake_args#*--configfile }

fi

configfile=${configfile%% *}
configfile=${configfile%% --configfile}
sptree=$(cat $configfile | shyaml get-value species_tree)
job_name=$(cat $configfile | shyaml get-value jobname)

#run SCORPiOs iteratively
for i in $(seq $iteration $max_iter); do

  #if it is not the first iteration, print names of corrected subtrees to file tmp_corrected_prev_iter
  if (( $i!=1 ))
    then
      cat "SCORPiOs_"${job_name}/Corrections/Accepted_Trees*$((i-1)) > .tmp_corrected_prev_iter_${job_name}
  else
    touch .tmp_corrected_prev_iter_${job_name}
  fi

  #run SCORPiOs if first iteration or if number of corrections in previous iter > min_correction
  if (( $i==1 )) || [[ $(wc -l <.tmp_corrected_prev_iter_${job_name}) -gt $min_correction ]]
    then
      echo "----------------"
      echo " Iteration: $i"
      echo "----------------"
      echo "Iteration: $i" >&2
      snakemake $snake_args $snake_config_args current_iter=$i --use-conda
      j=$i
  else
    break
  fi

done

#remove temp
rm .tmp_corrected_prev_iter_${job_name}

#If output exists (i.e no raised errors above) write .nhx correction tags and exit
if [ -f "SCORPiOs_${job_name}/SCORPiOs_output_${j}.nhx" ]; then

  echo "Termination after $j correction iterations ">&2
  echo "Browsing the corrected forests of each iteration to write final .nhx correction tags ">&2

  input="SCORPiOs_${job_name}/SCORPiOs_output_%d.nhx"
  output="SCORPiOs_${job_name}/SCORPiOs_output_${j}_with_tags.nhx"

  python -m scripts.trees.iteration_nhx_tags -o $output -i $j -c $input --internal -sp $sptree
  # python -m scripts.trees.iteration_nhx_tags -o $output -i $j -c $input
fi
