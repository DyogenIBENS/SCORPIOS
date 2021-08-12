import glob
import os
PWD = os.getcwd()

ctree_aore = "recombination_teleosts_small_v2/aore_ctrees"
ctree_lore = "recombination_teleosts_small_v2/lore_ctrees"
subtree_lore = "recombination_teleosts_small_v2/lore_trees"
subtree_aore = "recombination_teleosts_small_v2/aore_trees"
subtree_ml = "recombination_teleosts_small_v2/ml_trees"
OUT = "recombination_teleosts_small_v2/lktest"
subalis = "recombination_teleosts_small_v2/subalis_for_au"


FOREST = "../SCORPiOs/SCORPiOs_genofishv3_small/SCORPiOs_output_5_with_tags.nhx"
ALIS = "../SCORPiOs/data/genofish_v3_small/GenomicusV3_ali_small.fa.gz"
SPTREE = "../SCORPiOs/data/genofish_v3_small/GenomicusV3_sptree_small.nwk"

config["anc"] = config.get("anc", "Neopterygii")
config["outgr"] = "Lepisosteus.oculatus"
config["sp_to_keep"] = "-set Lepisosteus.oculatus Amia.calva Arapaima.gigas Danio.rerio Scleropages.formosus Paramormyrops.kingsleyae Oryzias.latipes Gasterosteus.aculeatus"

rule all:
    input: ".touch_autests"

# checkpoint subalis_loretrees_aoretrees:
#     input:
#         forest = FOREST,
#         sptree = SPTREE,
#         ali = ALIS
#     output:
#         alis = directory(f"{subalis}"),
#         trees_lore = directory(f"{ctree_lore}"),
#         trees_aore = directory(f"{ctree_aore}"),
#     params: anc = config["anc"], outgr = config["outgr"], args_restrict = config.get("sp_to_keep", '')
#     shell:
#         "python src/constrained_aore_lore_topologies.py -t {input.forest} -s {input.sptree} "
#         "--anc {params.anc} -o {output.alis} -ol {output.trees_lore} -oa {output.trees_aore} "
#         "-sp {params.outgr} -a {input.ali} {params.args_restrict}"


rule check_ali:
    input: f"{subalis}/{{tree}}.fa"
    output: f"{subalis}/{{tree}}.reduced.fa"
    conda: "../SCORPiOs/envs/raxml.yaml"
    shell: 
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.tree} -m GTRGAMMA "
        "-s {input} -w {PWD}/{subalis}/; "
        "if [ ! -s {output} ]; then cp {input} {output}; fi; "
        "rm {subalis}/RAxML_info.{wildcards.tree}"

rule lore_tree:
    input: ali = f"{subalis}/{{tree}}.reduced.fa", ctree = f"{ctree_lore}/{{tree}}.nh"
    output: subtree = f"{subtree_lore}/{{tree}}.nh"
    conda: "../SCORPiOs/envs/raxml.yaml"
    shell:
        "raxmlHPC -g {input.ctree} -n {wildcards.tree}_lore -m GTRGAMMA -p 1234 "
        "-s {input.ali} -w {PWD}/{subalis}/; "
        "mv {subalis}/RAxML_bestTree.{wildcards.tree}_lore {output}"

rule aore_tree:
    input: ali = f"{subalis}/{{tree}}.reduced.fa", ctree = f"{ctree_aore}/{{tree}}.nh"
    output: subtree = f"{subtree_aore}/{{tree}}.nh"
    conda: "../SCORPiOs/envs/raxml.yaml"
    shell:
        "raxmlHPC -g {input.ctree} -n {wildcards.tree}_aore -m GTRGAMMA -p 1234 "
        "-s {input.ali} -w {PWD}/{subalis}/; "
        "mv {subalis}/RAxML_bestTree.{wildcards.tree}_aore {output}"

rule ml_tree:
    input: ali = f"{subalis}/{{tree}}.reduced.fa"
    output: subtree = f"{subtree_ml}/{{tree}}.nh"
    conda: "../SCORPiOs/envs/raxml.yaml"
    shell:
        "raxmlHPC -p 1234 -n {wildcards.tree}_ml -m GTRGAMMA -s {input.ali} -w {PWD}/{subalis}/; "
        "mv {subalis}/RAxML_bestTree.{wildcards.tree}_ml {output}"

rule lk_test:
    input:
        ali = f"{subalis}/{{tree}}.reduced.fa", ml = f"{subtree_ml}/{{tree}}.nh",
        aore= f"{subtree_aore}/{{tree}}.nh", lore = f"{subtree_lore}/{{tree}}.nh"
    output: OUT+"/Res_{tree}.txt"
    threads: 1
    conda: "../SCORPiOs/envs/raxml.yaml"
    shell:
        "bash src/prototype_au_test3.sh {wildcards.tree} {input.ml} {input.ali} {input.aore} "
        "{input.lore} {output} || touch {output};"


def get_result(wildcards):
    #check that checkpoint has been executed
    # co = checkpoints.subalis_loretrees_aoretrees.get(**wildcards).output[0]
    subtrees, = glob_wildcards(subalis+"/{tree}.fa")
    subtrees = [i for i in subtrees if "reduced" not in i]
    out = expand(OUT+"/Res_{tree}.txt", tree=subtrees)
    return out


rule dummy:
    input: get_result
    output: touch(".touch_autests")
