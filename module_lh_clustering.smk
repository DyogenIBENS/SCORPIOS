import os

from ete3 import Tree

PWD = os.getcwd()

ALIS = SCORPIOS_CONFIG["alis"]

OUTFOLDER = f"SCORPiOs-LH_{JNAME}/clustering"

SP, GENES = list(config["dup_genome"].items())[0]
SORTBY = config.get("sort_by", "names")

arg_subset = config.get("sp_to_keep", '')
if arg_subset:
    arg_subset = "-set "+ arg_subset

arg_balance = config.get("balance_ohno", False)
if arg_balance:
    arg_balance = "--balance_ohno"

# rule Target:
#     input:
#         expand(f"{JOBNAME}/hmm_clusters_on_genome_{{i}}clusters_{{j}}states.svg", i=range(2, 5),
#                j=range(2, 5))

checkpoint extract_subtrees_subalis:
    input:
        forest = SCORPIOS_CORRTREES,
        sptree = SPTREE,
        ali = ALIS
    output:
        trees = directory(f"{OUTFOLDER}/subtrees/"),
        alis = directory(f"{OUTFOLDER}/subalis/")
    params: anc = LORE_WGD, outgr = f"{LORE_OUTGR}", arg_subset = arg_subset, arg_balance = arg_balance
    shell:
        "python -m scripts.lore_hunter.sample_subtrees_for_clust -t {input.forest} -s {input.sptree} "
        "--anc {params.anc} -oa {output.alis} -ot {output.trees} -sp {params.outgr} "
        "-a {input.ali} {params.arg_subset} {params.arg_balance}"

rule check_ali:
    input: f"{OUTFOLDER}/subalis/{{subtrees}}.fa"
    output: f"{OUTFOLDER}/subalis/{{subtrees}}.reduced.fa"
    # conda: "envs/raxml.yaml"
    shell: 
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.subtrees} -m GTRGAMMA "
        "-s {input} -w {PWD}/{OUTFOLDER}/subalis/; "
        "if [ ! -s {output} ]; then cp {input} {output}; fi; "
        "rm {OUTFOLDER}/subalis/RAxML_info.{wildcards.subtrees}"

rule prepare_tree:
    input: f"{OUTFOLDER}/subtrees/{{subtrees}}.nhx"
    output: f"{OUTFOLDER}/subtrees/{{subtrees}}.reduced.nhx"
    run:
        t = Tree(input[0], format=1)
        t.write(outfile=output[0], format=9)

rule brlength:
    input:
        tree = f"{OUTFOLDER}/subtrees/{{subtrees}}.reduced.nhx",
        ali = f"{OUTFOLDER}/subalis/{{subtrees}}.reduced.fa"
    output: f"{OUTFOLDER}/subtrees/{{subtrees}}_raxml.nhx"
    # conda: "envs/raxml.yaml"
    shell:
        "raxmlHPC -f e -n {wildcards.subtrees} -m GTRGAMMA -s {input.ali} -t {input.tree} "
        "-w {PWD}/{OUTFOLDER}/subalis/; mv {OUTFOLDER}/subalis/RAxML_result.{wildcards.subtrees} {output}; "
        "rm {OUTFOLDER}/subalis/RAxML_*.{wildcards.subtrees}"

rule sp_tag:
    input: t = f"{OUTFOLDER}/subtrees/{{subtrees}}.nhx", b = f"{OUTFOLDER}/subtrees/{{subtrees}}_raxml.nhx"
    output: f"{OUTFOLDER}/subtrees/{{subtrees}}_final.nhx"
    run:
        t1 = Tree(input[0])
        t2 = Tree(input[1])
        for leaf in t1.get_leaves():
            sp = leaf.S
            t2_leaf = t2.get_leaves_by_name(name=leaf.name)[0]
            t2_leaf.S = sp

        outgr1, outgr2 = t1.get_children()
        outgr = outgr1
        if len(outgr1) > len(outgr2):
            outgr = outgr2

        if outgr1.is_leaf():
            t2.set_outgroup(outgr1.name)

        elif outgr2.is_leaf():
            t2.set_outgroup(outgr2.name)

        else:
            anc = t2.get_common_ancestor({i.name for i in outgr1.get_leaves()})
            if anc == t2:
                anc = t2.get_common_ancestor({i.name for i in outgr2.get_leaves()})
            t2.set_outgroup(anc)

        #reroot (raxml unroots my fixed topology)
        t2.write(outfile=output[0], format=1, features=['S'])


def get_all_subtrees(wildcards):
    #check that checkpoint has been executed
    co = checkpoints.extract_subtrees_subalis.get(**wildcards).output[0]
    subtrees, = glob_wildcards(co+"/{subtree}.nhx")
    subtrees = [i for i in subtrees if "reduced" not in i and "raxml" not in i and '_final' not in i]
    out = expand(OUTFOLDER + "/subtrees/{subtrees}_final.nhx", subtrees=subtrees)
    return out


rule dummy:
    input: get_all_subtrees
    output: touch(f"{OUTFOLDER}/.touch")

rule tree_distances:
    input: t = OUTFOLDER+".touch", incons = SUMMARY
    output: f"{OUTFOLDER}/dist_mat.csv"
    conda: "envs/tree_dist.yaml"
    params: trees = OUTFOLDER + "/subtrees/"
    threads: 50
    shell:
        "python src/trees_distances.py -t {params.trees} -i {input.incons} -nc {threads} -o {output}; "
        "cat {OUTFOLDER}/dist_mat*.csv > {output}"

#TODO add medoid representation of clusters and of each consistent/inconsistent group (top 10 medoids?)
rule treeclust:
    input: f"{OUTFOLDER}/dist_mat.csv"
    output:
        expand(OUTFOLDER+"/clust_{i}.tsv", i=range(2, 5)),
        expand(OUTFOLDER+"/mds_{i}.svg", i=range(2, 5)),
        expand(OUTFOLDER+"/mat_{i}.svg", i=range(2, 5))
    conda: "envs/treecl.yaml"
    params: clust = OUTFOLDER+"/clust.tsv", mds = OUTFOLDER+"/mds.svg", mat = OUTFOLDER+"/mat.svg"
    threads: 50
    shell:
        "python src/tree_clust.py -d {input} -o {params.clust} -om {params.mds} -od {params.mat} -nc {threads}"

rule prepare_clusters_for_genome_plot:
    input: c = f"{OUTFOLDER}/clust_{{i}}.tsv", t = OUTFOLDER+".touch"
    output: f"{OUTFOLDER}/clust_{{i}}_for_plot.tsv"
    params: trees = OUTFOLDER + "/subtrees/", sp = "Arapaima.gigas" #TODO: put in config
    shell:
        "python src/clustering_res_to_dupsp.py -c {input.c} -t {params.trees} -s {params.sp} -o {output}"

# rule consensus_topologies:
#     input: clusters = f"{OUTFOLDER}/clust_{{i}}.tsv", trees = get_all_subtrees
#     output: lambda wildcards: expand(OUTFOLDER+"/consensus_tree_{j}.nhx", j=range(1, wildcards.i))
#     params: odir = f"SCORPiOs-LH_{OUTFOLDER}/"
#     shell: "python scripts.lore_hunter.consensus_topologies -c {input.clusters} -t {input.trees} -o {params.odir}"

rule plot_clusters:
    input:
        fam = f"{OUTFOLDER}/clust_{{i}}_for_plot.tsv",
        genes = GENES
        
    output: fig = f"{OUTFOLDER}/clusters_on_genome_{{i}}.png", data = f"{OUTFOLDER}/clusters_on_genome_{{i}}.pkl"
    params: sp = SP
    shell:
        "python ../paralogy_map/src/plot_paralogy_map.py -c {input.fam} -g {input.genes} -o {output.fig} "
        "-s {params.sp} -sort {SORTBY} -f dyogen -t 'tree topology clusters (k={wildcards.i})' --default_palette --singlesp --save"


rule fit_hmm:
    input:
        clust = f"{OUTFOLDER}/clusters_on_genome_{{i}}.pkl"        
    output: fig = f"{OUTFOLDER}/hmm_clusters_on_genome_{{i}}clusters_{{j}}states.svg"
    params: sp = SP
    conda: "envs/hmm.yaml"
    shell:
        "python src/fit_simple_hmm.py -c {input.clust} -o {output.fig} -n {wildcards.j} --sp {params.sp}"