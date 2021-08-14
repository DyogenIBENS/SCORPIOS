from snakemake.io import load_configfile
import sys

#may need --scheduler greedy


#get scorpios config,
SCORPIOS_CONFIGFILE = config["scorpios_config"]

#SCORPiOs is an extension of SCORPiOs, and depends on SCORPiOs as a subworkflow.
subworkflow scorpios:
    workdir:
        '.'
    configfile:
        SCORPIOS_CONFIGFILE


#Set all output names FIXME (duplicate code from main scorpios)
def out_name(name, jobname, iteration, wcard_wgd='', wcard_outgr=''):
    """
    Generates output names with jobname directory prefix and iteration suffix.
    Also adds required wildcards.
    """
    if wcard_wgd:
        name+="_"+wcard_wgd

    if wcard_outgr:
        name+="_"+wcard_outgr

    name = "SCORPiOs_"+jobname+'/'+name+'_'+str(iteration)
    return name

#check lh mode
assert config.get("mode", "diagnostic").lower() in ["clustering", "likelihood_tests", "diagnostic"],\
       "Invalid `mode`, please check your config."

MODE = config["mode"].lower()

REC_BR = config.get("recompute_brln", False)
assert isinstance(REC_BR, bool)

REC_ALI = config.get("recompute_msa", False)
assert isinstance(REC_ALI, bool)

# check that configs are consistent
SCORPIOS_CONFIG = load_configfile(SCORPIOS_CONFIGFILE)
if len(SCORPIOS_CONFIG["WGDs"]) == 1:
    LORE_WGD = list(SCORPIOS_CONFIG["WGDs"].keys())[0]
    assert (not config.get("lore_wgd", "") or config.get("lore_wgd", "") == LORE_WGD)

else:
    assert config.get("lore_wgd", ""), "SCORPiOs was run to correct several WGDs, you should\
                                        specify in config the one you want to run LORe analyses on."
    LORE_WGD = config["lore_wgd"]

LORE_OUTGR = config.get("lore_outgr", "")
if not LORE_OUTGR:
    LORE_OUTGR = SCORPIOS_CONFIG["WGDs"][LORE_WGD].split(",")
    if len(LORE_OUTGR) > 1:
        sys.stderr.write(f"Selecting {LORE_OUTGR[0]} as outgroup.")
    LORE_OUTGR = LORE_OUTGR[0]

#Get SCORPiOs inputs and outputs
JNAME = SCORPIOS_CONFIG["jobname"]
ITER = config.get("iter", 0)
CTREES = scorpios(out_name("Trees/ctrees", JNAME, ITER))
Acc = scorpios(out_name("Corrections/Accepted_Trees", JNAME, ITER, LORE_WGD))
ORTHOTABLE = scorpios(out_name("Families/Homologs", JNAME, ITER, LORE_WGD, LORE_OUTGR))
SUMMARY = scorpios(out_name("Corrections/Trees_summary", JNAME, ITER, LORE_WGD))
RES = scorpios(out_name("Corrections/Res_polylk", JNAME, ITER))
SCORPIOS_CORRTREES = scorpios(out_name("SCORPiOs_output", JNAME, ITER)+'.nhx')

SPTREE = SCORPIOS_CONFIG["species_tree"]
if "trees" in SCORPIOS_CONFIG and ITER == 0:
    INPUT_FOREST = SCORPIOS_CONFIG["trees"]
elif "trees" not in SCORPIOS_CONFIG and ITER == 0:
    INPUT_FOREST = out_name("input_forest", JNAME, ITER)+'.nhx'
if ITER != 0:
    INPUT_FOREST = out_name("SCORPiOs_output", JNAME, ITER) + '.nhx'

print(Acc)
print(CTREES)
print(SUMMARY)

CTREES_DIR = scorpios(CTREES+"/"+LORE_WGD+"/")


### WORKFLOW

print(MODE)
if MODE.lower() == "diagnostic":
    rule Target:
        input:
            f"SCORPiOs-LH_{JNAME}/diagnostic/seq_synteny_conflicts_by_homeologs.svg",
            f"SCORPiOs-LH_{JNAME}/diagnostic/seq_synteny_conflicts_on_genome.svg"

elif MODE.lower() == "clustering":
    rule Target:
        input:
            expand("SCORPiOs-LH_"+JNAME+"/clustering/medoids_clustering_k_"+str(config.get("k", 3))+"/cluster_{i}_medoid_{n}.nhx",
                   i=range(0, config.get("k", 3)), n=range(1, config.get("n", 5)+1)),
            expand("SCORPiOs-LH_"+JNAME+"/clustering/medoids_incons/medoid_{n}.nhx",
                   n=range(1, config.get("n", 5)+1)),
            "SCORPiOs-LH_"+JNAME+"/clustering/inconsistent_trees_vs_clusters.txt",
            "SCORPiOs-LH_"+JNAME+"/clustering/clusters_k-"+str(config.get("k", 3))+"_on_genome.svg"

else:
    rule Target:
        input:
            "SCORPiOs-LH_"+JNAME+"/lktests/lore_aore_on_genome.svg"

rule check_scorpios_output_integrity:
    input: scorpios("SCORPiOs_"+JNAME+"/.cleanup_"+str(ITER))
    output: touch(f"SCORPiOs-LH_{JNAME}/integrity_checkpoint.out")
    run: 
        ctrees_a, = glob_wildcards(CTREES+"/"+LORE_WGD+"/C_{ctrees}.nh")
        ctrees_b, = glob_wildcards(RES+"/"+LORE_WGD+"/Res_{ctrees}.txt")
        sys.stderr.write('Checking SCORPiOs output integrity...\n')
        if not set(ctrees_a) and set(ctrees_a) == set(ctrees_b):
            print("Please re-run scorpios, output of the checkpoint rule appears to be incomplete.")
            sys.exit(1)

def get_ctrees(wildcards, restrict=None):
    Ctrees, = glob_wildcards(CTREES+"/"+LORE_WGD+"/C_{ctrees}.nh")
    Ctrees = [i.split('/')[-1] for i in Ctrees]
    if restrict:
        Ctrees = [i for i in Ctrees if restrict in i]

    out = expand(CTREES+"/"+LORE_WGD+"/C_{ctrees}.nh", ctrees=Ctrees)
    # out = expand('test_lore/{ctrees}_test.txt', ctrees=Ctrees)
    return out


#include the 3 SCORPiOs LORe Hunter modules
include: "module_lh_diagnostic.smk"
include: "module_lh_clustering.smk"
include: "module_lh_au-test.smk"
