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

N = config.get("n", 5)
K = config.get("k", 3)


print(OUTFOLDER)
#TODO: linting

#FIXME: handle outgroup better to have correspondance with the Accepted SCORPiOs file
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


def get_all_subtrees(wildcards, OUTFOLDER=OUTFOLDER):
    #check that checkpoint has been executed
    co = checkpoints.extract_subtrees_subalis.get(**wildcards).output[0]
    subtrees, = glob_wildcards(co+"/{subtree}.nhx")
    subtrees = [i for i in subtrees if "reduced" not in i and "raxml" not in i and '_final' not in i]
    out = expand(OUTFOLDER + "/subtrees/{subtrees}_final.nhx", subtrees=subtrees)
    return out


rule dummy:
    input: get_all_subtrees
    output: touch(f"{OUTFOLDER}/.touch")


rule prepare_summary:
    input: ctreedir = scorpios(CTREES_DIR), summary = scorpios(SUMMARY), acc = scorpios(Acc)
    output: fam = f"{OUTFOLDER}/trees_summary.txt"
    shell:
        "python -m scripts.lore_hunter.write_ancgenes_treeclust -a {input.acc} -t {input.ctreedir} "
        "--summary_only -c {input.summary} -o {output.fam}"

rule tree_distances:
    input: t = OUTFOLDER+"/.touch", incons = f"{OUTFOLDER}/trees_summary.txt"
    output: temp(f"{OUTFOLDER}/distance_matrix.csv")
    conda: "envs/tree_dist.yaml"
    params: trees = OUTFOLDER + "/subtrees/", odir = OUTFOLDER
    threads: 50 #max_cores in config
    shell:
        "python -m scripts.lore_hunter.trees_distances -t {params.trees} -i {input.incons} -nc {threads} -o {params.odir}; "
        "cat {OUTFOLDER}/dist_mat*.csv > {output}; rm {OUTFOLDER}/dist_mat*.csv"

rule treeclust:
    input: f"{OUTFOLDER}/distance_matrix.csv"
    output:
        clust = f"{OUTFOLDER}/clust_{K}.tsv",
        mds = f"{OUTFOLDER}/mds_{K}.svg",
        mat = f"{OUTFOLDER}/mat_{K}.svg",
        omat = f"{OUTFOLDER}/distance_matrix_converted.csv"
    conda: "envs/treecl.yaml"
    params: k = K
    shell:
        "python -m scripts.lore_hunter.tree_clust -d {input} -o {output.clust} -om {output.mds} "
        "-od {output.mat} -k {params.k}"

rule medoids_clust:
    input: dmat = f"{OUTFOLDER}/distance_matrix_converted.csv", clust = f"{OUTFOLDER}/clust_{K}.tsv"
    output: expand(OUTFOLDER+'/medoids_clustering_k_'+str(K)+"/cluster_{i}_medoid_{n}.nhx", i=range(0, K), n=range(1, N+1))
    params: n = N, odir = OUTFOLDER+'/medoids_clustering_k_'+str(K), trees = OUTFOLDER + "/subtrees/{}_final.nhx"
    shell: "python -m scripts.lore_hunter.medoid_topologies -d {input.dmat} -c {input.clust} -n {params.n} -t {params.trees} -o {params.odir}"

rule medoids_incons:
    input: dmat = f"{OUTFOLDER}/distance_matrix_converted.csv", incons = f"{OUTFOLDER}/trees_summary.txt"
    output: expand(OUTFOLDER+'/medoids_incons/medoid_{n}.nhx', n=range(1, N+1))
    params: n = N, trees = OUTFOLDER + "/subtrees/{}_final.nhx", odir = OUTFOLDER+'/medoids_incons',
    shell: "python -m scripts.lore_hunter.medoid_topologies -d {input.dmat} -c {input.incons} -n {params.n} -t {params.trees} -o {params.odir} -r 'Inconsistent'"

rule compare_incons_clusters:
    input: clust = f"{OUTFOLDER}/clust_{K}.tsv", incons = f"{OUTFOLDER}/trees_summary.txt"
    output: OUTFOLDER+'/inconsistent_trees_vs_clusters.txt'
    shell: "python -m scripts.lore_hunter.compare_incons_clusters -c {input.clust} -i {input.incons} -o {output}"

rule clusters_full_summary:
    input: treedir = f"{OUTFOLDER}/subtrees/", clusters = f"{OUTFOLDER}/clust_{K}.tsv"
    output: f"{OUTFOLDER}/clustering_summary_k-{K}.tsv"
    shell: "python -m scripts.lore_hunter.write_ancgenes_treeclust -t {input.treedir} "
        "-c {input.clusters} -o {output}"

rule prepare_clusters_for_rideogram:
    input: c = f"{OUTFOLDER}/clustering_summary_k-{K}.tsv", genes = GENES
    output: karyo = f"{OUTFOLDER}/karyo_ide.txt",
            feat = f"{OUTFOLDER}/clust_k-{K}_ide.txt"        
    shell:
        "python -m scripts.lore_hunter.make_rideograms_inputs -i {input.c} -g {input.genes} "
        "-k {output.karyo} -o {output.feat} -f dyogen"

rule plot_clusters_on_genome:
    input:
        karyo = f"{OUTFOLDER}/karyo_ide.txt",
        feat = f"{OUTFOLDER}/clust_k-{K}_ide.txt"
    output: temp(f"{OUTFOLDER}/clusters_k-{K}_on_genome_tmp.svg")
    params: k=K
    conda: 'envs/rideogram.yaml'
    shell:
        "Rscript scripts/lore_hunter/plot_genome.R -k {input.karyo} -f {input.feat} -o {output} -c {params.k}"

rule rm_legend_rideogram:
    input: f"{OUTFOLDER}/clusters_k-{K}_on_genome_tmp.svg"
    output: temp(f"{OUTFOLDER}/clusters_k-{K}_on_genome_tmp2.svg")
    shell: "sed 's/Low.*//g' {input} | sed 's/\\(.*\\)\\<text.*/\\1\\/svg\\>/' > {output}"

def make_title(default_title, input_name):
    if 'hmm' in input_name:
        return default_title + "(hmm fit)"
    else:
        return default_title

rule add_legend_and_title_rideogram:
    input: f = f"{OUTFOLDER}/clusters_k-{K}_on_genome_tmp2.svg"
    output: f"{OUTFOLDER}/clusters_k-{K}_on_genome.svg"
    params:
        k=K,
        title=f'Tree topologies clusters on {SP} chromosomes',
        labels=expand("'cluster {i}'", i=range(0, K))
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lore_hunter.fix_rideogram -i {input} -o {output} -c {params.k} "
        "-t '{params.title}' -l {params.labels}"

# if config.get("fit_hmm", False):
#     #TODO: lk test to find best fitting HMM
#     rule fit_hmm:
#         input:
#             clust = f"{OUTFOLDER}/clust_k-{K}_ide.txt"       
#         output: fig = f"{OUTFOLDER}/clustering/clusters_k-{K}_hmm_fit_ide.txt"
#         params:
#             sp = SP
#         conda: "envs/hmm.yaml"
#         shell:
#             "python src/fit_simple_hmm.py -c {input.clust} -o {output.fig} -n {wildcards.j} --sp {params.sp}"


#     rule plot_hmm_on_genome:
#         input:
#             karyo = f"{OUTFOLDER}/karyo_ide.txt",
#             feat = f"{OUTFOLDER}/clust_k-{K}_hmm_fit_ide.txt"
#         output: temp(f"{OUTFOLDER}/clusters_k-{K}_hmm_fit_on_genome_tmp.svg")
#         params: k=k
#         conda: 'envs/rideogram.yaml'
#         shell:
#             "Rscript scripts/lore_hunter/plot_genome.R -k {input.karyo} -f {input.feat} -o {output} -c {params.k}"