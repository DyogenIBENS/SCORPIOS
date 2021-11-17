
"""
SCORPiOs LORelEi is an extension to the main SCORPiOs pipeline, it analyzes sequence-synteny
conflicts in gene tree.
"""

from snakemake.io import load_configfile
import sys


# Get SCORPiOs config
SCORPIOS_CONFIGFILE = config["scorpios_config"]

# Update main SCOPRiOs current_iter if necessary
ITERATION = config.get("iter", 0)

if ITERATION != 0:
    SCORPIOS_CONFIGFILE_COPY = '.'.join(SCORPIOS_CONFIGFILE.split(".")[:-1]) + '.copy.yaml'
    with open(SCORPIOS_CONFIGFILE, 'r') as infile, open(SCORPIOS_CONFIGFILE_COPY, 'w') as outfile:
        iter_updated = False
        for line in infile:
            if 'current_iter' in line:
                outfile.write(f"current_iter: {ITERATION}\n")
                iter_updated = True
            else:
                outfile.write(line)
        if not iter_updated:
            outfile.write(f"current_iter: {ITERATION}\n")

    SCORPIOS_CONFIGFILE = SCORPIOS_CONFIGFILE_COPY

# SCORPiOs LORelEi is an extension of SCORPiOs, and depends on SCORPiOs as a subworkflow.
subworkflow scorpios:
    workdir:
        '.'
    configfile:
        SCORPIOS_CONFIGFILE


# Set all output names
def out_name(name, JOBNAME_S, iteration, wcard_wgd='', wcard_outgr=''):
    """
    Generates output names with JOBNAME_S directory prefix and iteration suffix.
    Also adds required wildcards.
    """
    if wcard_wgd:
        name+="_"+wcard_wgd

    if wcard_outgr:
        name+="_"+wcard_outgr

    name = "SCORPiOs_"+JOBNAME_S+'/'+name+'_'+str(iteration)
    return name

# Check LORelEi mode
assert config.get("mode", "diagnostic").lower() in ["likelihood_tests", "diagnostic"],\
       "Invalid `mode`, please check your config."

MODE = config["mode"].lower()

# Check that configs are consistent
SCORPIOS_CONFIG = load_configfile(SCORPIOS_CONFIGFILE)

if len(SCORPIOS_CONFIG["WGDs"]) == 1:
    LORE_WGD = list(SCORPIOS_CONFIG["WGDs"].keys())[0]
    assert (not config.get("lore_wgd", "") or config.get("lore_wgd", "") == LORE_WGD),\
           "Invalid `lore_wgd`, please check your config."

else:
    assert config.get("lore_wgd", ""), "SCORPiOs was run to correct several WGDs, you should\
                                        specify in config the one you want to run LORe analyses on."
    LORE_WGD = config["lore_wgd"]


LORE_OUTGRS = SCORPIOS_CONFIG["WGDs"][LORE_WGD]
LORE_OUTGR = LORE_OUTGRS.split(',')[0]

if MODE == "diagnostic":
    if len(LORE_OUTGRS.split(',')) > 1:
        assert "use_anc" in config["pre_dup_proxy"] or "use_outgr" in config["pre_dup_proxy"]
        if "use_outgr" in config["pre_dup_proxy"]:
            LORE_OUTGR = config["pre_dup_proxy"]["use_outgr"]
    else:
        config["pre_dup_proxy"] = config.get("pre_dup_proxy", {})
        config["pre_dup_proxy"]["use_outgr"] = config.get("use_outgr", LORE_OUTGR)
    

# Get SCORPiOs outputs
JOBNAME_S = SCORPIOS_CONFIG["jobname"]

CONSTREES = scorpios(out_name("Trees/ctrees", JOBNAME_S, ITERATION))
ACCEPTED = scorpios(out_name("Corrections/Accepted_Trees", JOBNAME_S, ITERATION, LORE_WGD))
ORTHOTABLE = scorpios(out_name("Families/Homologs", JOBNAME_S, ITERATION, LORE_WGD, LORE_OUTGR))

ORTHOTABLES = []
for OUTGROUP in LORE_OUTGRS.split(','):
    ORTHOTABLES.append(scorpios(out_name("Families/Homologs", JOBNAME_S, ITERATION, LORE_WGD, OUTGROUP)))
SUMMARY = scorpios(out_name("Corrections/Trees_summary", JOBNAME_S, ITERATION, LORE_WGD))
RES = scorpios(out_name("Corrections/Res_polylk", JOBNAME_S, ITERATION))
SCORPIOS_CORRTREES = scorpios(out_name("SCORPiOs_output", JOBNAME_S, ITERATION)+'.nhx')
COMBIN = out_name("Graphs/outcombin", JOBNAME_S, ITERATION, LORE_WGD)

SPTREE = SCORPIOS_CONFIG["species_tree"]

CTREES_DIR = scorpios(CONSTREES+"/"+LORE_WGD+"/")

JOBNAME_L = JOBNAME_S
if "jname" in config:
    JOBNAME_L += '_' + config["jname"]


# LORELEI CONFIG (ALL modes)

SP = config["dup_genome"]
GENES = SCORPIOS_CONFIG["genes"] % SP


# Set LORelEi WORKFLOW Targets

if MODE.lower() == "diagnostic":
    rule Target:
        input:
            f"SCORPiOs-LORelEi_{JOBNAME_L}/diagnostic/seq_synteny_conflicts_by_homeologs.svg",
            f"SCORPiOs-LORelEi_{JOBNAME_L}/diagnostic/seq_synteny_conflicts_on_genome.svg"

else:
    rule Target:
        input:
            "SCORPiOs-LORelEi_"+JOBNAME_L+"/lktests/lore_aore_on_genome.svg"

rule check_scorpios_output_integrity:
    """
    Explicitly verifies that intermediary outputs from SCORPiOs are complete.
    """
    input: scorpios("SCORPiOs_"+JOBNAME_S+"/.cleanup_"+str(ITERATION))
    output: touch(f"SCORPiOs-LORelEi_{JOBNAME_L}/integrity_checkpoint.out")
    run: 
        ctrees_a, = glob_wildcards(CONSTREES+"/"+LORE_WGD+"/C_{ctrees}.nh")
        ctrees_b, = glob_wildcards(RES+"/"+LORE_WGD+"/Res_{ctrees}.txt")
        sys.stdout.write(f"SCORPiOs LORelEi: {MODE} mode\n")
        sys.stderr.write('Checking SCORPiOs output integrity...\n')
        if not set(ctrees_a) or set(ctrees_a) != set(ctrees_b):
            print("Please re-run SCORPiOs, output of the checkpoint rule appears to be incomplete.")
            sys.exit(1)


#include SCORPiOs LORelEi modules
include: "module_lorelei_diagnostic.smk"
include: "module_lorelei_lktests.smk"
