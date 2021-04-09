PATTERN = config["pattern"] #TODO --> remove this somehow or have a default like 'chr*'
SP, GENES = list(config["dup_genome"].items())[0]
SORTBY = config.get("sort_by", "names")
assert "use_anc" in config.get("pre_dup_proxy", "") or "use_outgr" in config.get("pre_dup_proxy", "")

if "use_outgr" in config["pre_dup_proxy"]:

    REF = config["pre_dup_proxy"]["use_outgr"]

    rule prepare_homeologs_outgr:
        input: fam = ORTHOTABLE, summary = SUMMARY, check = f"SCORPiOs-LH_{JNAME}/integrity_checkpoint.out"
        output: incons = f"SCORPiOs-LH_{JNAME}/conflicts", all_trees = f"SCORPiOs-LH_{JNAME}/trees"
        shell:
            "grep {PATTERN} {input.fam} | cut -f 1 | uniq -c > {output.all_trees}; "
            "grep {PATTERN} {input.fam} | cut -f 1,3 > SCORPiOs-LH_{JNAME}/tmp_genes; "
            "grep Incons {input.summary} | cut -f 1 > SCORPiOs-LH_{JNAME}/tmp_incons; "
            "grep -f SCORPiOs-LH_{JNAME}/tmp_incons SCORPiOs-LH_{JNAME}/tmp_genes | cut -f 1 | uniq -c > {output.incons}"

else:

    REF_FILE = config["pre_dup_proxy"]["use_anc"]
    REF = "Pre-duplication"

    rule prepare_homeologs_anc: #TODO: fam could take several orthothables (ORTHOTABLES.split())
        input: fam = ORTHOTABLE, summary = SUMMARY, pm = REF_FILE, check = f"SCORPiOs-LH_{JNAME}/integrity_checkpoint.out"
        output: incons = f"SCORPiOs-LH_{JNAME}/conflicts", alltrees = f"SCORPiOs-LH_{JNAME}/trees"
        shell:
            "python -m scripts.lore_hunter.homeologs_pairs_from_paralogymap -i {input.fam} -p {input.pm} "
            "-s {input.summary} -oi {output.incons} -oa {output.alltrees}"

rule plot_homeologs:
    input: incons = f"SCORPiOs-LH_{JNAME}/conflicts", all_trees = f"SCORPiOs-LH_{JNAME}/trees"
    output: f"SCORPiOs-LH_{JNAME}/seq_synteny_conflicts_by_homeologs.svg"
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lore_hunter.homeologs_tree_conflicts -i {input.incons} -g {input.all_trees} -o {output} "
        " --refname '{REF}'"

rule prepare_genome_plot:
    input: ctreedir = CTREES_DIR, summary = SUMMARY
    output: fam = f"SCORPiOs-LH_{JNAME}/inconsistent_families.tsv", pal = f"SCORPiOs-LH_{JNAME}/palette.tsv"
    shell:
        "grep Incons {input.summary} | grep -v multigenic > SCORPiOs-LH_{JNAME}/incons; echo -e 'Inconsistent\tred' > {output.pal}; "
        "python scripts/lore_hunter/write_ancgenes_treeclust.py -t {input.ctreedir} -c SCORPiOs-LH_{JNAME}/incons -o {output.fam}" #RIdeogram + bundle adjacent


#TODO: if possible use RIdeogram
#TODO: if possible highlight high-densty regions
#TODO: bundle adjacent genes together to decrease image size
rule plot_conflicts_on_genome:
    input:
        fam = f"SCORPiOs-LH_{JNAME}/inconsistent_families.tsv",
        pal = f"SCORPiOs-LH_{JNAME}/palette.tsv",
        genes = GENES
    output: f"SCORPiOs-LH_{JNAME}/seq_synteny_conflicts_on_genome.svg"
    params: sp = SP
    conda: "envs/plots.yaml"
    shell:
        "python -m scripts.lore_hunter.plot_genome -c {input.fam} -g {input.genes} -pf {input.pal} "
        "-o {output} -s {params.sp} -sort {SORTBY} -f dyogen -t 'sequence-synteny conflicts' --save"


rule get_conflicts_by_ancestors:
    input: ctrees = get_ctrees("ENSLOC"), ori_forest = INPUT_FOREST 
    output: f"SCORPiOs-LH_{JNAME}/conflicts_by_anc_by_sp.csv"
    shell:
        "touch {output}"
        # "python -m scripts.trees.inconsistent_trees.py --by_anc --no_ctrees --by_sp"

#TODO: get_conflicts_by_ancestor
# --> mettre à jour inconsistent_trees
# --> faire en sorte de pouvoir tester en enlevant les sp descendant de plusieurs ancêtres +  des sp_ind
# --> il faudra prbabalement faire des copies des objets arbres de ete3 sinon je vais casser le script et casser scorpios par la même occaz :/

rule plot_conflicts_by_ancestors:
    input: f"SCORPiOs-LH_{JNAME}/conflicts_by_anc_by_sp.csv"
    output: f"SCORPiOs-LH_{JNAME}/seq_synteny_conflicts_by_ancestors_and_sp.svg"
    shell: "touch {output}"