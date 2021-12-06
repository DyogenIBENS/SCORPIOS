import glob
import os
import shutil


PWD = os.getcwd()
ALIS = SCORPIOS_CONFIG["alis"]
OUTFOLDER = f"SCORPiOs-LORelEi_{JOBNAME_L}/lktests"

arg_subset = config.get("sp_to_keep", '')
if arg_subset:
    arg_subset = "-set "+ arg_subset

LORE_CLASSES = config.get("lore_groups", {"LORE": 'default'})
arg_groups = ''
if LORE_CLASSES["LORE"] != 'default':
    arg_groups = '-gr '+ LORE_CLASSES["LORE"]
LABELS = ' '.join(["AORe"] + list(LORE_CLASSES.keys()))

RAXML_SEED = config.get("raxml_seed", 1234)

MIN_SIZE = config.get("min_size", None)

ARG_MIN = ""
if MIN_SIZE is not None:
    ARG_MIN = f"--min_size {MIN_SIZE}"

checkpoint subalis_loretrees_aoretrees:
    """
    Extract WGD gene families and gene alignments from SCORPiOs corrected forest and builds
    corresponding constrained AORe and LORe gene tree topologies.
    """
    input:
        forest = SCORPIOS_CORRTREES,
        sptree = SPTREE,
        ali = ALIS,
        ctreedir = CTREES_DIR
    output:
        alis = directory(f"{OUTFOLDER}/subalis/"),
        trees_lore = directory(f"{OUTFOLDER}/ctree_lore/"),
        trees_aore = directory(f"{OUTFOLDER}/ctree_aore/"),
    params: anc = LORE_WGD, outgr = LORE_OUTGRS.replace(',', ' '), arg_subset = arg_subset, gr = arg_groups
    shell:
        "python -m scripts.lorelei.constrained_aore_lore_topologies -t {input.forest} -c {input.ctreedir} -s {input.sptree} "
        "--anc {params.anc} -o {output.alis} -ol {output.trees_lore} -oa {output.trees_aore} "
        "-sp {params.outgr} -a {input.ali} {params.arg_subset} {params.gr}"


rule check_ali_lktests:
    """
    Checks the input alignment with RAxML.
    """
    input: f"{OUTFOLDER}/subalis/{{tree}}.fa"
    output: temp(f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa")
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree} || true && "
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.tree} -m GTRGAMMA "
        "-s {input} -w {PWD}/{OUTFOLDER}/subalis/ >&2 && "
        "if [ ! -s {OUTFOLDER}/subalis/{wildcards.tree}.fa.reduced ]; "
        "then cp {input} {output}; else mv {OUTFOLDER}/subalis/{wildcards.tree}.fa.reduced {output};fi"


rule aore_lore_tree:
    """
    Builds the LORe and AORe (depending on the `class` wildcard) tree with RAxML.
    Logs are saved in {class}_trees/RAxML_info.{tree}.
    """
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ctree = f"{OUTFOLDER}/ctree_{{class}}/{{tree}}.nh"
    output: tree = f"{OUTFOLDER}/{{class}}_trees/{{tree}}.nh",
            log = f"{OUTFOLDER}/{{class}}_trees/RAxML_info.{{tree}}",
            tmp_log = temp(f"{OUTFOLDER}/subalis/RAxML_log.{{tree}}_{{class}}"),
            tmp_tree = temp(f"{OUTFOLDER}/subalis/RAxML_result.{{tree}}_{{class}}"),
    params: raxml_seed = RAXML_SEED
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_{wildcards.class} || true &&"
        "raxmlHPC -g {input.ctree} -n {wildcards.tree}_{wildcards.class} -m GTRGAMMA -p {params.raxml_seed} "
        "-s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/ >&2 && "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_{wildcards.class} {output.tree} && "
        "mv {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_{wildcards.class} {output.log}"


rule ml_tree:
    """
    Builds the Maximum Likelihood (ML) tree with RAxML. 
    Logs are saved in ml_trees/RAxML_info.{tree}.
    """
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa"
    output: tree = f"{OUTFOLDER}/ml_trees/{{tree}}.nh",
            log = f"{OUTFOLDER}/ml_trees/RAxML_info.{{tree}}",
            tmp_tree = temp(f"{OUTFOLDER}/subalis/RAxML_result.{{tree}}_ml"),
    params: raxml_seed = RAXML_SEED
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_ml || true &&"
        "raxmlHPC -p {params.raxml_seed} -n {wildcards.tree}_ml -m GTRGAMMA -s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/ && "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_ml {output.tree} >&2 && "
        "mv {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_ml {output.log}"


rule lk_test:
    """
    Runs the likelihood AU test to compare likelihoods of the AORe, LORe and ML trees. 
    """
    input:
        ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ml = f"{OUTFOLDER}/ml_trees/{{tree}}.nh",
        aore = f"{OUTFOLDER}/aore_trees/{{tree}}.nh", lore = f"{OUTFOLDER}/lore_trees/{{tree}}.nh"
    output: f"{OUTFOLDER}/lktests/Res_{{tree}}.txt"
    shell:
        "bash scripts/prototype_au_test3.sh {wildcards.tree} {input.ml} {input.ali} {input.aore} "
        "{input.lore} {output} && rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_lktest"


def get_result(wildcards):

    # Check that checkpoint has been executed
    co = checkpoints.subalis_loretrees_aoretrees.get(**wildcards).output[0]

    # Ouputs from the checkpoint
    subtrees, = glob_wildcards(OUTFOLDER+"/subalis/{tree}.fa")
    subtrees = [i for i in subtrees if ".reduced" not in i]

    # Expand expected outputs 
    out = expand(OUTFOLDER+"/lktests/Res_{tree}.txt", tree=subtrees)

    return out

rule list_lktest:
    """
    Lists all of the CONSEL Likelihood-tests result files (files to parse in later steps). 
    """
    input: get_result
    output: outf = temp(OUTFOLDER+"/file_list.txt")
    run:
        with open(output.outf, 'w') as fw1:
            for f in input:
                fw1.write(f+'\n')

rule make_summary:
    """
    Writes a short summary of likelihood-tests results confronting the LORe and AORe hypotheses.
    Lists the family IDs (outgroup gene name) of AORe and LORe gene trees.
    """
    input: OUTFOLDER+"/file_list.txt"
    output: OUTFOLDER+"/lore_aore_summary_au_all.txt"
    shell:
        "python -m scripts.trees.parse_au_test -i {input} -o {output} --lore -w {LORE_WGD}"


rule clean_subalis_folder:
    input: OUTFOLDER+"/lore_aore_summary_au_all.txt"
    output: touch(OUTFOLDER+"/.clean_alis")
    params: clean = config.get("clean_alis", True)
    run:
        if params.clean:
            try:
                shutil.rmtree(f"{OUTFOLDER}/subalis/")
            except OSError as e:
                print(f"Error: {e.filename} - {e.strerror}.")

rule lore_aore_full_summary:
    """
    Writes a full summary of likelihood-tests results confronting the LORe and AORe hypotheses.
    Lists all genes in the AORe and LORe gene trees.
    """
    input: clusters = OUTFOLDER+"/lore_aore_summary.txt", clean = OUTFOLDER+"/.clean_alis"
    output: f"{OUTFOLDER}/lore_aore_summary.tsv"
    params: treedir = f"{OUTFOLDER}/ml_trees/" #FIXME: may break stuff for reruns (explicit input?)
    shell: "python -m scripts.lorelei.write_ancgenes_treeclass -t {params.treedir} "
           "-c {input.clusters} -o {output} -r 'lore rejected' 'aore rejected'"


rule prepare_lore_aore_for_rideogram:
    """
    Prepares input files for the RIdeogram karyotype plot.
    """
    input: c = f"{OUTFOLDER}/lore_aore_summary.tsv", genes = GENES
    output: karyo = f"{OUTFOLDER}/karyo_ide.txt",
            feat = f"{OUTFOLDER}/lore_aore_ide.txt"
    params: arg_min = ARG_MIN
    shell:
        "python -m scripts.lorelei.make_rideograms_inputs -i {input.c} -g {input.genes} "
        "-k {output.karyo} -o {output.feat} -f dyogen {params.arg_min}"


rule plot_lore_aore_on_genome:
    """
    Uses RIdeograms to plot aore and lore gene families on the karyotype of a duplicated genome.
    """
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/lore_aore_ide.txt"
    output: temp(f"{OUTFOLDER}/lore_aore_on_genome_tmp.svg")
    params: nb_classes = len(LABELS.split())
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lorelei/plot_genome.R -k {input.karyo} -f {input.feat} -o {output} -c {params.nb_classes}"


rule rm_legend:
    """
    Removes automatically generated legend from RIdeogram (because it assumes a continuous variable)
    """
    input: f"{OUTFOLDER}/lore_aore_on_genome_tmp.svg"
    output: temp(f"{OUTFOLDER}/lore_aore_on_genome_tmp2.svg")
    shell: "sed 's/Low.*//g' {input} > {output} && echo '</text></svg>' >> {output}"


rule add_legend_and_title:
    """
    Adds a correct legend and a title for the RIdeogram plot.
    """
    input: f"{OUTFOLDER}/lore_aore_on_genome_tmp2.svg"
    output: f"{OUTFOLDER}/lore_aore_on_genome.svg"
    params: sp = SP, labels = LABELS, nb_classes = len(LABELS.split())
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lorelei.fix_rideogram -i {input} -o {output} -c {params.nb_classes} "
        "-t 'AORe and LORe topologies on {params.sp} chromosomes' -l {params.labels}"