
"""
SCORPiOs snakemake module to identify orthologous gene communtities in orthology graphs and derive
topological constraints on corresponding gene trees.
"""
SPECTRAL = config.get('spectral', '')
if SPECTRAL and SPECTRAL.lower() not in ['n', 'no', 'false']:
    SPECTRAL = "--spectral"
else:
    SPECTRAL = ""

SUBSET_TO_CORRECT = config.get('subset_to_correct', '') #should be a file, one gene per line (one gene per tree to correct), if provided in config
if SUBSET_TO_CORRECT:
    SUBSET_TO_CORRECT = '--correct_only '+SUBSET_TO_CORRECT

rule cut_orthology_graphs:
    """
    Loads orthology graphs and uses community detection algorithms to identify the 2 orthogroups.
    """
    input: Sorted_SyntenyOrthoPred+'.gz'
    output: a=GraphsOrthogroups, b=Summary
    threads: config["ncores"]
    params: Summary = Summary.replace("{{outgr}}", "{{wildcards.outgr}}")\
                                     .replace("{{wgd}}","{{wildcards.wgd}}")
    conda: "envs/graphs.yaml"
    shell:"""
    python -m scripts.graphs.orthogroups -i {input} -o {output.a} -n {threads} -s {params.Summary}\
    -ignSg {config[ignoreSingleGeneCom]} -wgd {wildcards.wgd},{wildcards.outgr} {SPECTRAL}
    """


def expand_graphs_outgr(wildcards):
    """
    Gets a list of graph cut results files (one for each outgroup)
    """
    return expand(GraphsOrthogroups.replace("{wgd}", "{{wgd}}"),\
                  outgr=config['WGDs'][wildcards.wgd].split(','))


def expand_graphcut_summary(wildcards):
    """
    Gets a list of graph cut summary files (one for each outgroup)
    """
    return expand(Summary.replace("{wgd}", "{{wgd}}"),\
                  outgr=config['WGDs'][wildcards.wgd].split(','))


def expand_orthotables_outgr(wildcards):
    """
    Gets a list of orthology tables (one for each outgroup)
    """
    return expand(OrthoTable.replace("{wgd}", "{{wgd}}"),\
                  outgr=config['WGDs'][wildcards.wgd].split(','))


checkpoint gene_trees_to_correct:
    """
    Converts orthogroups into topological constraints on the gene tree.
    Save all subtrees, subalignment and constrained tree where contraint is absent in input trees.
    """
    input: graph_cuts=expand_graphs_outgr, graph_cuts_summaries=expand_graphcut_summary,
           orthotables=expand_orthotables_outgr, alis = config['alis']

    output: ctrees = directory(CTREES+"/{wgd}/"), subalis = directory(SUBALIS+"/{wgd}/"),
            subtrees = directory(SUBTREES+"/{wgd}/"), tsum = TREES_SUMMARY

    params: graphs = lambda wildcards, input: ",".join(list(input.graph_cuts)),
            gsum = lambda wildcards, input: ",".join(list(input.graph_cuts_summaries)),
            otable = lambda wildcards, input: ",".join(list(input.orthotables)),
            outgroups = lambda wildcards: config['WGDs'][wildcards.wgd],
            outcombin = outcombin.replace("{{wgd}}", "{{wildcards.wgd}}"),
            args_subset_to_correct = SUBSET_TO_CORRECT

    shell:"""
    python -m scripts.trees.inconsistent_trees -n {params.outgroups} -i {params.graphs}\
    -f {params.otable} -t {input_trees} -a {config[alis]} -oc {output.ctrees} -oa {output.subalis}\
    -ot {output.subtrees} -gs {params.gsum} -s {output.tsum} -wgd {wildcards.wgd}\
    -fcombin {params.outcombin} {params.args_subset_to_correct}
    """
