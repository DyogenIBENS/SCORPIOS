#!/usr/bin/env bash

name=$1
enstree=$2
alifile=$3
cortree=$4
output=$5


workingDir=$(pwd)

alidir=${alifile%/*}
ensdir=${enstree%/*}
cordir=${cortree%/*}


#check alignment and remove undetermined columns (if any)
raxmlHPC -f c --print-identical-sequences -n ${name} -m GTRGAMMA --HKY85 -s "${alifile}" -w "${workingDir}/${alidir}/" >&2
if [ -s "${alifile}.reduced" ]; then
	alifile="${alifile}.reduced"
fi
rm "${alidir}/RAxML_info.${name}"


#if not already computed for original tree, compute sites likelihood under HKY model for both trees
if [ ! -s "${alidir}/${name}_a.lk" ]; then

	#cat the two tree to a single file and remove nhx tags
	cat "${cortree}" "${enstree}" | sed -e 's/\[[^][]*\]//g' > "${cordir}/trees_${name}.nh"

	#compute site lk
	raxmlHPC -f G -n ${name} -m GTRGAMMA --HKY85 -s "${alifile}" -z "${cordir}/trees_${name}.nh" -w "${workingDir}/${alidir}/" >&2
	rm "${cordir}/trees_${name}.nh"

	#rename output
	mv "${alidir}/RAxML_perSiteLLs.${name}" "${alidir}/${name}_a.lk"
	cp "${alidir}/${name}_a.lk" "${alidir}/tmp_${name}.lk"

#compute sites likelihood under HKY model for new solution tree only (original already computed)
else
	#remove nhx tags
	sed -e 's/\[[^][]*\]//g' "${cortree}" > "${cordir}/tmp_${name}.nh"

	#compute site lk
	raxmlHPC -f g -n ${name} -m GTRGAMMA --HKY85 -s "${alifile}" -z "${cordir}/tmp_${name}.nh" -w "${workingDir}/${alidir}/" >&2
	rm "${cordir}/tmp_${name}.nh"

	#extract original tree site log-lk and put it together in one file
	sed '1h;1d;'3'x' <(sed 2'!d;q' "${alidir}/RAxML_perSiteLLs.${name}") "${alidir}/${name}_a.lk" > "${alidir}/${name}_b.lk"
	cp "${alidir}/${name}_b.lk" "${alidir}/tmp_${name}.lk"

fi


#workaround since consel decides to trim filenames containing '.' (looks like extension split issue)
namenew="${name%.*}"

#test if difference in likelihood is signifiant with the AU-Test using consel
makermt --puzzle "${alidir}/tmp_${name}.lk" >&2
consel "${alidir}/tmp_${name}" >&2
catpv "${alidir}/tmp_$namenew.pv" > "$output"

## CLEAN ALL TEMP ##
#clean all consel temp
rm "${alidir}/tmp_${name}.lk"
rm "${alidir}/tmp_${name}.rmt"
rm "${alidir}/tmp_${name}.vt"
rm "${alidir}/tmp_${name}"
rm "${alidir}/tmp_${namenew}.pv"
rm "${alidir}/tmp_${namenew}.ci"

#clean raxml temp (log, already saved above)
rm "${alidir}/RAxML_info.${name}"
if [ -s "${alifile}.reduced" ]; then
	rm "${alifile}.reduced"
fi