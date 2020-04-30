
"""
SCORPiOs snakemake module to build starting gene trees with TreeBeST.
"""

if ITER < 2: #allows to use --forceall in iterative mode without this rule executing at iter 2

    rule build_input_trees:
        """
        Build input trees with TreeBeST best.
        """
        input: a = config["alis"], sp = config["species_tree"], map = config["genes_sp_mapping"]
        output: input_trees
        threads: config['ncores']
        conda: "envs/treebest_raxml_consel.yaml"
        shell:"""
        python -m scripts.trees.build_treebest_trees -a {input.a} -s {input.sp} -m {input.map}\
                                                     -nc {threads} -o {output}\
                                                     -tmp SCORPiOs_{config[jobname]}/
        """
