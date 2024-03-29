
"""
SCORPiOs snakemake module to predict gene orthologies and paralogies using synteny conservation
patterns, for all pairs of duplicated species.
"""

## parallelization scheme --> snakemake parallelization is inefficient with a large number of jobs
## if many duplicated species, I convert lots of small jobs into a smaller number of bigger jobs
args_pairwise_parallel = ""
BIG_RUN = False
if "parallel_scheme_large_job" in config and config["parallel_scheme_large_job"] == 'y':
    BIG_RUN = True
    args_pairwise_parallel = "-chr_list "+Chr


rule synteny_orthologies_paralogies:
    """
    Identifies post-WGD orthologs and paralogs in all pairs of duplicated species, using synteny.
    Executed in parallel for each chromsome in the outgroup and each pair of duplicated species.
    """
    input: a=OrthoTableF, d=TreesOrthologies+'/{wgd}', e=regions
    output: temp(Pairwise_SyntenyOrthoPred+"/{pairwise}_{chr_outgr}_{outgr}_{wgd}.txt")
    params: args_autho = args_autho.replace("{{outgr}}", "{{wildcards.outgr}}").replace("{{wgd}}",\
                                            "{{wildcards.wgd}}"),
            args_parallel = args_pairwise_parallel.replace("{{outgr}}", "{{wildcards.outgr}}")\
                                                  .replace("{{wgd}}","{{wildcards.wgd}}")
    shell:"""
    python -m scripts.synteny.pairwise_orthology_synteny -i {input.a} \
    -p {wildcards.pairwise} -chr {wildcards.chr_outgr} \
    -ortho {input.d} -o {output} -w {config[windowSize]} \
    -cutoff={config[cutoff]} {params.args_autho} {params.args_parallel}
    """


def all_output_species_pair_outgr_chr(wildcards):
    """
    Expands the chr_outgr wildcards and returns the list of expected outputs of rule
    synteny_orthologies_paralogies.
    """
    chr_outgr = []
    with open(checkpoints.outgroup_chromosomes.get(**wildcards).output[0]) as infile:

        chr_outgr += [line.strip() for line in infile]

    if BIG_RUN:

        chr_outgr = ["all_chrom"]

    ALL_DUPLICATED_SPECIES = spt.get_species(config["species_tree"], wildcards.wgd,
                                             ','.join(config["WGDs"].keys()), lowcov)
    ALL_DUPLICATED_SPECIES = list(itertools.combinations(ALL_DUPLICATED_SPECIES, 2))
    ALL_SPECIES_PAIRS = [('_').join((i, j)) for (i, j) in ALL_DUPLICATED_SPECIES]

    return expand(Pairwise_SyntenyOrthoPred+"/{pairwise}_{chr_outgr}_{{outgr}}_{{wgd}}.txt",
                  pairwise=ALL_SPECIES_PAIRS, chr_outgr=chr_outgr)


rule orthology_graphs:
    """
    Concatenates all pairwise outputs into a single file.
    Note: A simple 'cat' command could fail with error 'Argument list too long', if too many inputs.
    """
    input: all_output_species_pair_outgr_chr
    output: o=temp(SyntenyOrthoPred)
    run:
        with open(output.o, 'w') as outfile:
            for fname in input:
                if wildcards.outgr in fname:
                    with open(fname) as infile:
                        outfile.write(infile.read())


rule orthology_graphs_sort_gzip:
    """
    Sorts and compress the file with all concatenated orthology, for more efficient loading.
    """
    input: SyntenyOrthoPred
    output: Sorted_SyntenyOrthoPred + '.gz'
    params: out = Sorted_SyntenyOrthoPred.replace("{{outgr}}", "{{wildcards.outgr}}")\
                                                  .replace("{{wgd}}", "{{wildcards.wgd}}"),
            tmp = "SCORPiOs_"+config["jobname"]+'/'
    shell:"""
    sort -k 4 {input} -o {params.out} -T {params.tmp} && gzip {params.out};
    """
