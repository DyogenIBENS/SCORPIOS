
"""
SCORPiOs LORelEi module to analyze the genomic location of sequence-synteny conflicted families?
"""

OUTFOLDER = f"SCORPiOs-LORelEi_{JOBNAME_L}/diagnostic"


OUTGR_GENES = SCORPIOS_CONFIG["genes"] % LORE_OUTGR
POST_DUP = config.get("is_post_dup", False)
if POST_DUP:
    POST_DUP = '--post_dup'
else:
    POST_DUP = ''

COMBIN_ARG = ''
if len(LORE_OUTGRS.split(',')) > 1:
    COMBIN_ARG += '-c '+ COMBIN

if "use_outgr" in config["pre_dup_proxy"]:

    REF = config["pre_dup_proxy"]["use_outgr"]

    rule prepare_homeologs_outgr:
        input:
            fam = ORTHOTABLE,
            summary = SUMMARY,
            acc = ACCEPTED,
            check = f"SCORPiOs-LORelEi_{JOBNAME_L}/integrity_checkpoint.out"
        output:
            incons = f"{OUTFOLDER}/conflicts",
            alltrees = f"{OUTFOLDER}/trees",
        params:
            genes = OUTGR_GENES,
            fomt = SCORPIOS_CONFIG.get("genes_format", "bed"),
            combin_args = COMBIN_ARG
        shell:
            "python -m scripts.lorelei.homeologs_pairs_from_ancestor -i {input.fam} --is_outgroup "
            "-homeo {params.genes} -s {input.summary} -oi {output.incons} -oa {output.alltrees} "
            "-a {input.acc} {params.combin_args} -f {params.fomt}"

else:

    REF_FILE = config["pre_dup_proxy"]["use_anc"]
    REF = "Pre-duplication"

    rule prepare_homeologs_anc: 
        input:
            fam = ORTHOTABLES,
            summary = SUMMARY,
            acc = ACCEPTED,
            pm = REF_FILE,
            check = f"SCORPiOs-LORelEi_{JOBNAME_L}/integrity_checkpoint.out"
        output: incons = f"{OUTFOLDER}/conflicts", alltrees = f"{OUTFOLDER}/trees"
        params: postdup = POST_DUP
        shell:
            "python -m scripts.lorelei.homeologs_pairs_from_ancestor -i {input.fam} -a {input.acc} "
            "-homeo {input.pm} -s {input.summary} -oi {output.incons} -oa {output.alltrees} "
            "{params.postdup}"

rule plot_homeologs:
    input: incons = f"{OUTFOLDER}/conflicts", all_trees = f"{OUTFOLDER}/trees"
    output: f"{OUTFOLDER}/seq_synteny_conflicts_by_homeologs.svg"
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lorelei.homeologs_tree_conflicts -i {input.incons} -g {input.all_trees} "
        "-o {output} --refname '{REF}'"

rule prepare_genome_plot:
    input: ctreedir = CTREES_DIR, summary = SUMMARY, acc = ACCEPTED,
    output: fam = f"{OUTFOLDER}/inconsistent_families.tsv"
    shell:
        "python -m scripts.lorelei.write_ancgenes_treeclass -a {input.acc} -t {input.ctreedir} "
        "-c {input.summary} -o {output.fam} -r 'Inconsistent'"

rule prepare_input_rideogram:
    input:
        fam = f"{OUTFOLDER}/inconsistent_families.tsv",
        genes = GENES
    output:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/incons_ide.txt"
    params: sp = SP
    shell: "python -m scripts.lorelei.make_rideograms_inputs -i {input.fam} -g {input.genes} "
           "-k {output.karyo} -o {output.feat} -f dyogen"

rule plot_conflicts_on_genome:
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/incons_ide.txt"
    output: temp(f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp.svg")
    params: sp = SP
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lorelei/plot_genome.R -k {input.karyo} -f {input.feat} -o {output}"

rule remove_legend_rideogram:
    input: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp.svg"
    output: temp(f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp2.svg")
    shell: "sed 's/Low.*//g' {input} > {output} && echo '</text></svg>' >> {output}"

rule new_legend_and_title_rideogram:
    input: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome_tmp2.svg"
    output: f"{OUTFOLDER}/seq_synteny_conflicts_on_genome.svg"
    params: sp = SP
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lorelei.fix_rideogram -i {input} -o {output} -c 1 "
        "-t 'Sequence-synteny conflicts on {params.sp} chromosomes' -l 'inconsistent trees'"

