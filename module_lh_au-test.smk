import glob
import os
PWD = os.getcwd()


ALIS = SCORPIOS_CONFIG["alis"]

OUTFOLDER = f"SCORPiOs-LH_{JNAME}/lktests"

arg_subset = config.get("sp_to_keep", '')
if arg_subset:
    arg_subset = "-set "+ arg_subset

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
        "python -m scripts.lore_hunter.constrained_aore_lore_topologies -t {input.forest} -c {input.ctreedir} -s {input.sptree} "
        "--anc {params.anc} -o {output.alis} -ol {output.trees_lore} -oa {output.trees_aore} "
        "-sp {params.outgr} -a {input.ali} {params.arg_subset}"


rule check_ali_lktests:
    input: f"{OUTFOLDER}/subalis/{{tree}}.fa"
    output: f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa"
    shell: 
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.tree} -m GTRGAMMA "
        "-s {input} -w {PWD}/{OUTFOLDER}/subalis/; "
        "if [ ! -s {output} ]; then cp {input} {output}; fi; "
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.tree}"

rule aore_lore_tree:
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ctree = f"{OUTFOLDER}/ctree_{{class}}/{{tree}}.nh"
    output: subtree = f"{OUTFOLDER}/{{class}}_trees/{{tree}}.nh"
    params: raxml_seed = config.get("raxml_seed", 1234)
    shell:
        "raxmlHPC -g {input.ctree} -n {wildcards.tree}_{wildcards.class} -m GTRGAMMA -p {params.raxml_seed} "
        "-s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/; "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_{wildcards.class} {output}"


rule ml_tree:
    input: ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa"
    output: subtree = f"{OUTFOLDER}/ml_trees/{{tree}}.nh"
    params: raxml_seed = config.get("raxml_seed", 1234)
    shell:
        "raxmlHPC -p {params.raxml_seed} -n {wildcards.tree}_ml -m GTRGAMMA -s {input.ali} -w {PWD}/{OUTFOLDER}/subalis/; "
        "mv {OUTFOLDER}/subalis/RAxML_bestTree.{wildcards.tree}_ml {output}"

rule lk_test:
    input:
        ali = f"{OUTFOLDER}/subalis/{{tree}}.reduced.fa", ml = f"{OUTFOLDER}/ml_trees/{{tree}}.nh",
        aore = f"{OUTFOLDER}/aore_trees/{{tree}}.nh", lore = f"{OUTFOLDER}/lore_trees/{{tree}}.nh"
    output: f"{OUTFOLDER}/lktest/Res_{{tree}}.txt"
    shell:
        "bash scripts/prototype_au_test3.sh {wildcards.tree} {input.ml} {input.ali} {input.aore} "
        "{input.lore} {output} || touch {output};"


def get_result(wildcards):
    #check that checkpoint has been executed
    co = checkpoints.subalis_loretrees_aoretrees.get(**wildcards).output[0]
    subtrees, = glob_wildcards(OUTFOLDER+"/subalis/{tree}.fa")
    subtrees = [i for i in subtrees if "reduced" not in i]
    out = expand(OUTFOLDER+"/lktest/Res_{tree}.txt", tree=subtrees)
    return out


# rule get_res:
#     input: get_result
#     output: touch(".touch_autests")


rule list_lktest:
    input: get_result
    output: outf = temp(OUTFOLDER+"/file_list.txt")
    run:
        with open(output.outf,'w') as fw1:
            for f in input:
                fw1.write(f+'\n')


rule make_summary:
    input: OUTFOLDER+"/file_list.txt"
    output: OUTFOLDER+"/lore_aore_summary.txt"
    shell:
        "python -m scripts.trees.parse_au_test -i {input} -o {output} --lh -w {LORE_WGD}"


rule lore_aore_full_summary:
    input: treedir = f"{OUTFOLDER}/ml_trees/", clusters = OUTFOLDER+"/lore_aore_summary.txt"
    output: f"{OUTFOLDER}/lore_aore_summary_ancgenes.tsv"
    shell: "python -m scripts.lore_hunter.write_ancgenes_treeclust -t {input.treedir} "
           "-c {input.clusters} -o {output} -r 'lore rejected' 'aore rejected'"

rule prepare_lore_aore_for_rideogram:
    input: c = f"{OUTFOLDER}/lore_aore_summary_ancgenes.tsv", genes = GENES
    output: karyo = f"{OUTFOLDER}/karyo_ide.txt",
            feat = f"{OUTFOLDER}/lore_aore_ide.txt"        
    shell:
        "python -m scripts.lore_hunter.make_rideograms_inputs -i {input.c} -g {input.genes} "
        "-k {output.karyo} -o {output.feat} -f dyogen"

# #TODO: fix colors and legend
# #TODO: add a title with sp name
rule plot_lore_aore_on_genome:
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/lore_aore_ide.txt"
    output: f"{OUTFOLDER}/lore_aore_on_genome.svg"
    params: sp = SP
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lore_hunter/plot_genome.R -k {input.karyo} -f {input.feat} -o {output} -c 2"