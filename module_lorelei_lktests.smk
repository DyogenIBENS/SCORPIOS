import glob
import os
PWD = os.getcwd()


ALIS = SCORPIOS_CONFIG["alis"]

OUTFOLDER = f"SCORPiOs-LORelEi_{JNAME}/lktests"

arg_subset = config.get("sp_to_keep", '')
if arg_subset:
    arg_subset = "-set "+ arg_subset

LORE_CLASSES = config.get("lore_groups", ["LORE"])
LABELS = ' '.join(["AORe"] + list(LORE_CLASSES.keys()))

#TODO make it possible to test different speciation points (+default is first speciation point)
checkpoint subalis_loretrees_aoretrees:
    input:
        forest = SCORPIOS_CORRTREES,
        sptree = SPTREE,
        ali = ALIS,
        ctreedir = CTREES_DIR
    output:
        alis = directory(f"{OUTFOLDER}/subalis/"),
        trees_lore = directory(f"{OUTFOLDER}/ctree_lore/"),
        trees_aore = directory(f"{OUTFOLDER}/ctree_aore/"),
    params: anc = LORE_WGD, outgr = f"{LORE_OUTGR} Amia.calva", arg_subset = arg_subset
    shell:
        "python -m scripts.lorelei.constrained_aore_lore_topologies -t {input.forest} -c {input.ctreedir} -s {input.sptree} "
        "--anc {params.anc} -o {output.alis} -ol {output.trees_lore} -oa {output.trees_aore} "
        "-sp {params.outgr} -a {input.ali} {params.arg_subset}"


rule check_ali_lktests:
    input: f"{OUTFOLDER}/subalis/{{tree}}.fa"
    output: f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa"
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree} || true && "
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.tree} -m GTRGAMMA "
        "-s {input} -w {PWD}/{OUTFOLDER}/subalis/ && "
        "if [ ! -s {OUTFOLDER}/subalis/{wildcards.tree}.fa.reduced ]; "
        "then cp {input} {output}; else mv {OUTFOLDER}/subalis/{wildcards.tree}.fa.reduced {output};fi"

rule aore_lore_tree:
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ctree = f"{OUTFOLDER}/ctree_{{class}}/{{tree}}.nh"
    output: subtree = f"{OUTFOLDER}/{{class}}_trees/{{tree}}.nh"
    params: raxml_seed = config.get("raxml_seed", 1234)
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_{wildcards.class} || true &&"
        "raxmlHPC -g {input.ctree} -n {wildcards.tree}_{wildcards.class} -m GTRGAMMA -p {params.raxml_seed} "
        "-s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/ && "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_{wildcards.class} {output}"

rule ml_tree:
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa"
    output: subtree = f"{OUTFOLDER}/ml_trees/{{tree}}.nh"
    params: raxml_seed = config.get("raxml_seed", 1234)
    shell:
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_ml || true &&"
        "raxmlHPC -p {params.raxml_seed} -n {wildcards.tree}_ml -m GTRGAMMA -s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/ && "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_ml {output}"


rule lk_test:
    input:
        ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ml = f"{OUTFOLDER}/ml_trees/{{tree}}.nh",
        aore = f"{OUTFOLDER}/aore_trees/{{tree}}.nh", lore = f"{OUTFOLDER}/lore_trees/{{tree}}.nh"
    output: f"{OUTFOLDER}/lktests/Res_{{tree}}.txt"
    shell:
        "bash scripts/prototype_au_test3.sh {wildcards.tree} {input.ml} {input.ali} {input.aore} "
        "{input.lore} {output} && rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}_lktest"


def get_result(wildcards):
    #check that checkpoint has been executed
    co = checkpoints.subalis_loretrees_aoretrees.get(**wildcards).output[0]
    subtrees, = glob_wildcards(OUTFOLDER+"/subalis/{tree}.fa")
    subtrees = [i for i in subtrees if ".reduced" not in i]
    out = expand(OUTFOLDER+"/lktests/Res_{tree}.txt", tree=subtrees)
    return out

rule list_lktest:
    input: get_result
    output: outf = OUTFOLDER+"/file_list.txt"
    run:
        with open(output.outf, 'w') as fw1:
            for f in input:
                fw1.write(f+'\n')

rule make_summary:
    input: OUTFOLDER+"/file_list.txt"
    output: OUTFOLDER+"/lore_aore_summary.txt"
    shell:
        "python -m scripts.trees.parse_au_test -i {input} -o {output} --lore -w {LORE_WGD}"


rule lore_aore_full_summary:
    input: clusters = OUTFOLDER+"/lore_aore_summary.txt"
    output: f"{OUTFOLDER}/lore_aore_summary_ancgenes.tsv"
    params: treedir = f"{OUTFOLDER}/ml_trees/" #may break stuff for reruns
    shell: "python -m scripts.lorelei.write_ancgenes_treeclust -t {params.treedir} "
           "-c {input.clusters} -o {output} -r 'lore rejected' 'aore rejected'"

rule prepare_lore_aore_for_rideogram:
    input: c = f"{OUTFOLDER}/lore_aore_summary_ancgenes.tsv", genes = GENES
    output: karyo = f"{OUTFOLDER}/karyo_ide.txt",
            feat = f"{OUTFOLDER}/lore_aore_ide.txt"        
    shell:
        "python -m scripts.lorelei.make_rideograms_inputs -i {input.c} -g {input.genes} "
        "-k {output.karyo} -o {output.feat} -f dyogen"

rule plot_lore_aore_on_genome:
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/lore_aore_ide.txt"
    output: temp(f"{OUTFOLDER}/lore_aore_on_genome_tmp.svg")
    params: nb_classes = len(LABELS.split())
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lorelei/plot_genome.R -k {input.karyo} -f {input.feat} -o {output} -c {params.nb_classes}"

rule rm_legend:
    input: f"{OUTFOLDER}/lore_aore_on_genome_tmp.svg"
    output: temp(f"{OUTFOLDER}/lore_aore_on_genome_tmp2.svg")
    shell: "sed 's/Low.*//g' {input} | sed 's/\\(.*\\)\\<text.*/\\1\\/svg\\>/' > {output}"

rule add_legend_and_title:
    input: f"{OUTFOLDER}/lore_aore_on_genome_tmp2.svg"
    output: f"{OUTFOLDER}/lore_aore_on_genome.svg"
    params: sp = SP, labels = LABELS, nb_classes = len(LABELS.split())
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lorelei.fix_rideogram -i {input} -o {output} -c {params.nb_classes} "
        "-t 'AORe and LORe topologies on {params.sp} chromosomes' -l {params.labels}"