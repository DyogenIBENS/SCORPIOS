#!/usr/bin/env bash

name=$1
mltree=$2
alifile=$3
aoretree=$4
loretree=$5
output=$6

workingDir=$(pwd)

alidir=${alifile%/*}

raxmlHPC -f c --print-identical-sequences -n ${name}_lktest -m GTRGAMMA --HKY85 -s "${alifile}" -w "${workingDir}/${alidir}/" >&2
if [ -s "${alifile}.reduced" ]; then
    alifile="${alifile}.reduced"
fi
rm "${alidir}/RAxML_info.${name}_lktest"


#if not already computed for original tree, compute sites likelihood under HKY model for both trees
if [ ! -s "$output" ]; then

	#cat the trees to a single file and remove nhx tags and remove '()' around a single leaf
	cat "${mltree}" "${aoretree}" "${loretree}" | sed -e 's/\[[^][]*\]//g' -e 's/(\([^,]*\))/\1/g' > "${alidir}/trees_${name}.nh"

	#compute site lk
	raxmlHPC -f g -n ${name}_lktest -m GTRGAMMA --HKY85 -s "${alifile}" -z "${alidir}/trees_${name}.nh" -w "${workingDir}/${alidir}/" >&2
	rm "${alidir}/trees_${name}.nh"
    echo ${alifile}

	#rename output
	mv "${alidir}/RAxML_perSiteLLs.${name}_lktest" "${alidir}/${name}.lk"
fi

#workaround since consel decides to trim filenames containing '.' (looks like extension split issue)
namenew="${name%.*}"

#test if difference in likelihood is signifiant with the AU-Test using consel
makermt --puzzle "${alidir}/${name}.lk" >&2
consel "${alidir}/${name}" >&2
catpv "${alidir}/$namenew.pv" > "$output"


## CLEAN ALL TEMP ##
#clean all consel temp
rm "${alidir}/${name}.lk"
rm "${alidir}/${name}.rmt"
rm "${alidir}/${name}.vt"
rm "${alidir}/${name}"
rm "${alidir}/${namenew}.pv"
rm "${alidir}/${namenew}.ci"
