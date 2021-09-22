OUTGR_CHR_TO_PLOT = config.get("outgr_chroms_for_plot", "")
SP, GENES = list(config["dup_genome"].items())[0]
assert "use_anc" in config.get("pre_dup_proxy", "") or "use_outgr" in config.get("pre_dup_proxy", "")

OUTFOLDER = f"SCORPiOs-LH_{JNAME}/diagnostic"

if "use_outgr" in config["pre_dup_proxy"]:

    REF = config["pre_dup_proxy"]["use_outgr"]

    if OUTGR_CHR_TO_PLOT != '':
        rule restrict_to_outgr_chromosomes:
            input: fam = ORTHOTABLE
            output: f"{OUTFOLDER}/ORTHOTABLE_FOR_PLOTS"
            params: chr_file = OUTGR_CHR_TO_PLOT
            shell:
                "grep -Fwf {params.chr_file} {input.fam} > {output}"

    else:

        rule pass_outgr_chr_filter:
            input: fam = ORTHOTABLE
            output: touch(f"{OUTFOLDER}/ORTHOTABLE_FOR_PLOTS")

    rule prepare_homeologs_outgr:
        input:
            fam = f"{OUTFOLDER}/ORTHOTABLE_FOR_PLOTS",
            summary = SUMMARY,
            acc = Acc,
            check = f"SCORPiOs-LH_{JNAME}/integrity_checkpoint.out"
        output:
            incons = f"{OUTFOLDER}/conflicts",
            all_trees = f"{OUTFOLDER}/trees",
            tmp_acc = temp(f"{OUTFOLDER}/tmp_acc")
        shell:#FIXME: not sure this works with outfolder that way
            "cut -f 1 {input.fam} | uniq -c > {output.all_trees}; "
            "cut -f 1,3 {input.fam} > {OUTFOLDER}/tmp_genes; "
            "grep Incons {input.summary} | cut -f 1 > {OUTFOLDER}/tmp_incons; "
            "cut -f 1 {input.acc} > {OUTFOLDER}/tmp_acc; "
            "grep -v {OUTFOLDER}/tmp_acc {OUTFOLDER}/tmp_incons > {OUTFOLDER}/tmp_uncorr; "
            "grep -f {OUTFOLDER}/tmp_incons {OUTFOLDER}/tmp_genes | cut -f 1 | uniq -c > {output.incons}"

else:

    REF_FILE = config["pre_dup_proxy"]["use_anc"]
    REF = "Pre-duplication"

    #TODO: fam could take several orthothables (ORTHOTABLES.split()), if several outgroups
    #Check that the config arg with several outgroups can work with all extentions
    rule prepare_homeologs_anc: 
        input:
            fam = ORTHOTABLE,
            summary = SUMMARY,
            acc = Acc,
            pm = REF_FILE,
            check = f"SCORPiOs-LH_{JNAME}/integrity_checkpoint.out"
        output: incons = f"{OUTFOLDER}/conflicts", alltrees = f"{OUTFOLDER}/trees"
        shell:
            "python -m scripts.lore_hunter.homeologs_pairs_from_paralogymap -i {input.fam} -p {input.pm} "
            "-s {input.summary} -oi {output.incons} -oa {output.alltrees} -a {input.acc}"

rule plot_homeologs:
    input: incons = f"{OUTFOLDER}/conflicts", all_trees = f"{OUTFOLDER}/trees"
    output: f"{OUTFOLDER}/seq_synteny_conflicts_by_homeologs.svg"
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lore_hunter.homeologs_tree_conflicts -i {input.incons} -g {input.all_trees} -o {output} "
        " --refname '{REF}'"

rule prepare_genome_plot:
    input: ctreedir = CTREES_DIR, summary = SUMMARY, acc = Acc,
    output: fam = f"{OUTFOLDER}/inconsistent_families.tsv"
    shell:
        "python -m scripts.lore_hunter.write_ancgenes_treeclust -a {input.acc} -t {input.ctreedir} "
        "-c {input.summary} -o {output.fam} -r 'Inconsistent'"

rule prepare_input_rideogram:
    input:
        fam = f"{OUTFOLDER}/inconsistent_families.tsv",
        genes = GENES
    output:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/incons_ide.txt"
    params: sp = SP
    shell: "python -m scripts.lore_hunter.make_rideograms_inputs -i {input.fam} -g {input.genes} "
           "-k {output.karyo} -o {output.feat} -f dyogen"

rule plot_conflicts_on_genome:
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/incons_ide.txt",
    output: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp.svg"
    params: sp = SP
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lore_hunter/plot_genome.R -k {input.karyo} -f {input.feat} -o {output}"

rule remove_legend_rideogram:
    input: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp.svg"
    output: temp(f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp2.svg")
    shell: "sed 's/Low.*//g' {input} | sed 's/\\(.*\\)\\<text.*/\\1\\/svg\\>/' > {output}"

rule new_legend_and_title_rideogram:
    input: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp2.svg"
    output: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome.svg"
    params: sp = SP
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lore_hunter.fix_rideogram -i {input} -o {output} -c 1 "
        "-t 'Sequence-synteny conflicts on {params.sp} chromosomes' -l 'inconsistent trees'"


#Introducing the LH extension: quick usage and 
#Extension to scorpios and can be invoked with -s scorpios_lh.smk --> will first run scorpios and then the lh extension on the output (check it runs scorpios completely)
# --> we could make this one use rejected correction of iteration 0.
#For iterative run
#--> first run scorpios with the wrapper and then lh extension will ensure integrity of output and use it
#input data are same as scorpios but needs an additional configuration file 
#Detailed presentation of supported modes : diagnostic, au, treecl
#--> detail new config
#--> make a small example which includes 1 chr in outgroup salmon with mix AORe and LORe
#Add a documentation for the API :)
#Add a comment on how to change only one param in config (and update job name for scorpios classic)
#Fix my spacing issues