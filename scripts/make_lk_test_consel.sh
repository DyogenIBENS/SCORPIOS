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

#safeguard for fastdist potential errors : check that solution tree exist and is non-zero
if [ -s "$cortree" ] ; then

	#generate shorter gene ids for phylip conversion
	python -m scripts.trees.convert_ids -a $alifile -t $cortree $enstree

	#convert to phylip
	seqret "fasta::${alidir}/tmp_$name.fa" "phylip::${alidir}/tmp_$name.phy"

	#compute sites likelihood under HKY model for solution tree
	phyml -i "${alidir}/tmp_$name.phy" -u "${cordir}/tmp_${name}.nh" -o lr --print_site_lnl -b 0 --quiet  >&2

	#if not already computed, compute sites likelihood under HKY model for original tree
	if [ ! -s "${alidir}/tmp_${name}_c.phy_phyml_lk.txt" ]; then

		cp "${alidir}/tmp_$name.phy" "${alidir}/tmp_${name}_c.phy"
		phyml -i "${alidir}/tmp_${name}_c.phy" -u "${ensdir}/tmp_$name.nh" -o lr --print_site_lnl -b 0 --quiet >&2
		rm "${alidir}/tmp_${name}_c.phy"
		rm "${alidir}/tmp_${name}_c.phy_phyml_stats.txt"
		rm "${alidir}/tmp_${name}_c.phy_phyml_tree.txt"
	fi

	#concatenate sites likelihood into a single file
	cat "${alidir}/tmp_${name}.phy_phyml_lk.txt" "${alidir}/tmp_${name}_c.phy_phyml_lk.txt" > "${alidir}/$name.lk"

	#workaround since consel decides to trim filenames containing '.'
	namenew="${name%.*}"

	#test if difference in likelihood is signifiant with the AU-Test using consel
	makermt --phyml "${alidir}/$name.lk" >&2
	consel "${alidir}/$name" >&2
	catpv "${alidir}/$namenew.pv" > "$output"

	## CLEAN ALL TEMP ##
	#clean all consel temp
	rm "${alidir}/${name}.lk"
	rm "${alidir}/${name}.rmt"
	rm "${alidir}/${name}.vt"
	rm "${alidir}/${name}"
	rm "${alidir}/${namenew}.pv"
	rm "${alidir}/${namenew}.ci"

	#clean pyhml temp
	rm "${alidir}/tmp_${name}.phy_phyml_lk.txt"
	rm "${alidir}/tmp_${name}.phy_phyml_stats.txt"
	rm "${alidir}/tmp_${name}.phy_phyml_tree.txt"


	#clean gene id conversion temp
	rm "${alidir}/tmp_$name.fa"
	rm "${alidir}/tmp_$name.phy"
	rm "$cordir/tmp_${name}.nh"
	rm "$ensdir/tmp_$name.nh"

else
	>&2 echo "Tree building failed."
	touch $output
	touch "${alidir}/tmp_${name}_c.phy_phyml_lk.txt"
fi
