
"""
SCORPiOs snakemake module to build test and regraft synteny aware subtrees for inconsistent subtrees.
"""

RESUME = config.get("resume", "n")

RAXML = config.get("brlength_tool", '')
RAXML_ARG = ''
if RAXML:
    assert RAXML.lower() in ["raxml", "treebest phyml"], "Bad 'brlength_tool' parameter in config."
    if RAXML.lower() == "raxml":
        RAXML_ARG = '--raxml'


rule remove_anc_in_sptree:
    """
    Removes ancestor names from the species tree, to input it to profileNJ.
    """
    input: config["species_tree"]
    output: NO_ANC_TREE
    run:
        spt.remove_anc(config["species_tree"], NO_ANC_TREE)


rule distance_matrix:
    """
    For each subtree to correct, builds a distance matrix using treebest distmat.
    """
    input: SUBALIS+"/{wgd}/{ctrees}.fa"
    output: temp(tmp_matrix+"_{ctrees}.phy")
    shell:"""
    treebest distmat kimura {input} > {output}
    """


rule build_tree_polytomysolver:
    """
    For each subtree to correct, find a corrected subtree in agreement with the constraints,
    using profileNJ.
    Note: An empty file is generated if fastdist failed.
    """
    input: ctree = CTREES+"/{wgd}/C_{ctrees}.nh", matrix = tmp_matrix+"_{ctrees}.phy",
           sp_tree = NO_ANC_TREE, ali = SUBALIS+"/{wgd}/{ctrees}.fa"
    output: PolyS+"/{wgd}/{ctrees}.nh"
    conda: "envs/polytomysolver.yaml"
    shell:"""
    if [ -s {input.matrix} ]; then
        profileNJ -s {input.sp_tree} -g {input.ctree} --sep '_' -d {input.matrix} -o {output}\
                  --slimit 1 --spos postfix  >&2;
    else touch {output};
    fi;
    sed -i.bak 1d {output}; rm {output}.bak;
    """


rule test_polyS:
    """
    Test profileNJ solution against the original tree using a likelihood test.
    """
    input:  polys = PolyS+"/{wgd}/{ctrees}.nh", ali = SUBALIS+"/{wgd}/{ctrees}.fa",
            ori_tree = SUBTREES+"/{wgd}/{ctrees}.nh"
    output: OutPolylk+"/{wgd}/Res_{ctrees}.txt", SUBALIS+"/{wgd}/{ctrees}_a.lk"
    threads: 1
    shell:"""
    bash scripts/make_lk_test_consel.sh {wildcards.ctrees} {input.ori_tree} {input.ali}\
                                        {input.polys} {output}
    """


rule build_test_tree_treebest:
    """
    Builds and tests a treebest phyml solution. Will generate empty file if profileNJ solution was
    accepted by lk-tests.
    """
    input: polylk = OutPolylk+"/{wgd}/Res_{ctrees}.txt",
           hkyensembl = SUBALIS+"/{wgd}/{ctrees}_a.lk",
           ali = SUBALIS+"/{wgd}/{ctrees}.fa",
           ori_tree = SUBTREES+"/{wgd}/{ctrees}.nh",
           ctree = CTREES+"/{wgd}/C_{ctrees}.nh"
    output: lktest=OuttreeBlk+"/{wgd}/Res_{ctrees}.txt",
            tree = treeB+"/{wgd}/{ctrees}.nh"
    params: outgroups=lambda wildcards: config['WGDs'][wildcards.wgd]
    shell:"""
    bash scripts/correct_subtrees_treebest.sh {wildcards.ctrees} {input.ali} {input.ctree}\
    {input.ori_tree} {output.lktest} {input.polylk} {params.outgroups} {output.tree}\
    {config[species_tree]}
    """


def get_ctrees(wildcards):
    """
    Gets a list of all expected likelihood-test results output.
    """
    #check that checkpoint has been executed
    co = checkpoints.gene_trees_to_correct.get(**wildcards).output[0]

    #get all generated ctrees to expand wildcards
    Ctrees, = glob_wildcards(co+"C_{ctrees}.nh")
    poly = expand(OutPolylk+"/{wgd}/Res_{ctrees}.txt", ctrees=Ctrees, wgd=wildcards.wgd)
    treeb = expand(OuttreeBlk+"/{wgd}/Res_{ctrees}.txt", ctrees=Ctrees, wgd=wildcards.wgd)
    return poly+treeb


rule list_treeb_inputs_for_parsing:
    """
    Creates two files listing all likelihood-test result files.

    Note: Simply giving all inputs to a script could fail with error 'Argument list too long',
           if too many inputs.
    """
     input: get_ctrees
     output: p = temp("SCORPiOs_"+config['jobname']+'/outfile_poly_{wgd}_'+str(ITER)),
             t = temp("SCORPiOs_"+config['jobname']+'/outfile_treeb_{wgd}_'+str(ITER))
     run:
         with open(output.p,'w') as fw1, open(output.t,'w') as fw2:
             for f in input:
                if 'polylk' in f:
                     fw1.write(f+'\n')
                elif "treeB" in f:
                    fw2.write(f+'\n')


rule successfully_corrected_trees:
    """
    Parses likelihood test results and writes a file listing accepted corrections.
    """
     input:  polySlk="SCORPiOs_"+config['jobname']+'/outfile_poly_{wgd}_'+str(ITER),
             treeblk="SCORPiOs_"+config['jobname']+'/outfile_treeb_{wgd}_'+str(ITER)
     output: Acc
     shell:"""
     python -m scripts.trees.parse_au_test -i {input.polySlk} -o {output} -it {input.treeblk}\
                                           -w {wildcards.wgd} -p {PolyS},{treeB}
     """


rule correct_input_trees:
    """
    Re-graft corrected subtrees in the full gene trees forest.
    """
    input: acc = expand(Acc, wgd=config['WGDs'].keys())
    output: a = outTrees
    threads: config.get("limit_threads_for_branch_lengths", config["ncores"])
    params: wgds = ','.join(config["WGDs"].keys()),
            outgroups = '_'.join([config["WGDs"][i] for i in config["WGDs"].keys()]),
            input = lambda wildcards, input: ",".join(list(input.acc))
    shell:"""
    python -m scripts.trees.regraft_subtrees -t {input_trees} -a {config[alis]} \
    -acc {params.input} -o {output.a} -s {config[species_tree]} \
    -n {threads} -anc {params.wgds} -ogr {params.outgroups} \
    -tmp {outTmpTrees} -sa {config[save_tmp_trees]} {arg_brlength} -r {RESUME} {RAXML_ARG}
    """

if config['save_subtrees_lktest'] == 'n':
    rule clean_interm:
        """
        Dirty rule to remove some intermediary results files, if specified in config. Currently,
        (due to checkpoints ?), it seems the temp() flag is not working correctly (re-triggering of
        rules, early removal...).
        """
        input: outTrees
        output: touch("SCORPiOs_"+config['jobname']+"/.int_cleanup_"+str(ITER))
        params: path = "SCORPiOs_"+config["jobname"]+"/Trees"
        shell:"""
        rm -r {params.path} || true; rm -r {PolyS} || true; rm -r {treeB} || true;\
        rm -r {OutPolylk} || true; rm -r {OuttreeBlk} || true
        """
else:
    rule pass_clean_interm:
        """
        Pass intermediary files removal.
        """
        input: outTrees
        output: "SCORPiOs_"+config['jobname']+"/.int_cleanup_"+str(ITER)
        shell:"""
        touch {output}
        """

rule clean_up:
    """
    Dirty rule to clean up temp files at the end of the workflow. Currently, (due to checkpoints ?),
    it seems the temp() flag is not working correctly (re-triggering of rules, early removal...).
    """
    input: "SCORPiOs_"+config['jobname']+"/.int_cleanup_"+str(ITER)
    output: touch("SCORPiOs_"+config['jobname']+"/.cleanup_"+str(ITER))
    params: path = "SCORPiOs_"+config["jobname"]+"/Families"
    shell:"""
    rm -r {SUBALIS} || true; rm -r {SUBTREES} || true; rm {NO_ANC_TREE} || true;\
    rm -r {TreesOrthologies} || true; rm {params.path}/Chr* || true;\
    rm {params.path}/HomologsStrict* || true;\
    rm {params.path}/HomologsFilter* || true;
    rm {params.path}/UNCERTAIN* || true;
    rm {params.path}/orthologs* || true;
    rm {params.path}/paralogs* || true;
    rm -r {Pairwise_SyntenyOrthoPred} || true
    """
