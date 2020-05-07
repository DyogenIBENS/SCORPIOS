import os
import glob
import itertools
from scripts.trees import speciestree as spt
from scripts.synteny import filter_regions

"""
 SCORPiOs corrects orthologies and paralogy relationships in gene trees, so that they are consistent
 with a known WGD event. To do so, the method takes advantage of synteny conservation patterns.
"""

## Set all output names
def out_name(name, jobname, iteration, wcard_wgd=False, wcard_outgr=False):
    """
    Generates output names with jobname directory prefix and iteration suffix.
    Also adds required wildcards.
    """
    if wcard_wgd:
        name+='_{wgd}'

    if wcard_outgr:
        name+='_{outgr}'

    name = "SCORPiOs_"+jobname+'/'+name+'_'+str(iteration)
    return name

JNAME = config["jobname"]
ITER = config['current_iter']

if "trees" in config and os.path.isfile(config["trees"]):
    input_trees = config["trees"]
else:
    input_trees = out_name("input_forest", JNAME, ITER)+'.nhx'

OrthoTableStrict = out_name("Families/OrthoTableStrict", JNAME, ITER, True, True)
Chr =  out_name("Families/Chr", JNAME, ITER, True, True)
OrthoTable = out_name("Families/OrthoTable", JNAME, ITER, True, True)
OrthoTableF = out_name("Families/OrthoTableFilter", JNAME, ITER, True, True)
TreesOrthologies = out_name("TreesOrthologies", JNAME, ITER)
SyntenyOrthoPred = out_name("SyntenyOrthoPred", JNAME, ITER, True, True)
Sorted_SyntenyOrthoPred = out_name("Synteny/Sorted_SyntenyOrthoPred", JNAME, ITER, True, True)
GraphsOrthogroups = out_name("Graphs/GraphsOrthogroups", JNAME, ITER, True, True)
Summary = out_name("Graphs/Summary", JNAME, ITER, True, True)
outcombin = out_name("Graphs/outcombin", JNAME, ITER, True)
Acc = out_name("Corrections/Accepted_Trees", JNAME, ITER, True)
MULTIGENIC = out_name("Trees/Multigenic", JNAME, ITER, True)
TREES_SUMMARY = out_name("Corrections/Trees_summary", JNAME, ITER, True)
SUBTREES = out_name("Trees/subtrees", JNAME, ITER)
SUBALIS = out_name("Trees/subalis", JNAME, ITER)
CTREES = out_name("Trees/ctrees", JNAME, ITER)
Pairwise_SyntenyOrthoPred = out_name("Pairwise_SyntenyOrthoPred", JNAME, ITER)
PolyS = out_name("Corrections/PolyS", JNAME, ITER)
OutPolylk = out_name("Corrections/Res_polylk", JNAME, ITER)
outTrees = out_name("SCORPiOs_output", JNAME, ITER)+'.nhx'
treeB= out_name("Corrections/TreeB", JNAME, ITER)
OuttreeBlk=out_name("Corrections/Res_treeBlk", JNAME, ITER)
outTmpTrees=out_name("Corrections/tmp_whole_trees", JNAME, ITER)
UNCERTAIN = out_name("Families/UNCERTAIN", JNAME, ITER, True, True)
regions = out_name("Families/tmp_iter_updated_regions", JNAME, ITER, True, True)
NO_ANC_TREE = out_name("sptree_no_anc", JNAME, ITER)+'.nwk'
tmp_matrix = out_name("fastdist_mat", JNAME, ITER, True)
fam_no_graph = out_name("Families/Summary_fam_no_graph", JNAME, ITER, True, True)
orthologs = out_name("Families/orthologs", JNAME, ITER, True, True)
paralogs = out_name("Families/paralogs", JNAME, ITER, True, True)
Threshold = out_name("Families/threshold", JNAME, ITER, True, True)

# get correct input names and params in iterative mode
args_autho = ''
OrthoTable_prev = ''
Acc_prev = ''
incombin = ''
if int(ITER) > 1:
    Acc_prev = out_name("Corrections/Accepted_Trees", JNAME, int(ITER)-1, True, False)
    OrthoTable_prev = out_name("Families/OrthoTable", JNAME, int(ITER)-1, True, True)
    input_trees = out_name("SCORPiOs_output", JNAME, int(ITER)-1)+'.nhx'
    args_autho = '-filter '+regions
    incombin = out_name("Graphs/outcombin", JNAME, int(ITER)-1, True)


arg_brlength = '-br '+str(config['brlength'])
#if in iterative mode we force re-computation of branch-lengths
# if int(ITER) > 0:
#     arg_brlength = '-br y'

## Set parameters from config
if "genes_sp_mapping" not in config:
    config["genes_sp_mapping"] = ""

# Set genes file format if dyogen format is specified (otherwise bed is assumed)
if "genes_format" not in config:
    config["genes_format"] = "bed"

#all WGDs in trees
anc_other = ''
anc_arg = '-ow '+','.join(config['WGDs'].keys())

#"lowcoverage" species, place them in families but don't use them in the graphs
lowcov = ''
lowcov_arg = ''
if "lowcov_sp" in config:
    with open(config["lowcov_sp"], 'r') as f:
        lowcov = (',').join([sp.strip() for sp in f])
    lowcov_arg = '-l '+lowcov

wildcard_constraints:
    pairwise="[A-Za-z.0-9]+_[A-Za-z.0-9]+" #no underscore in sp names, as also required in ensembl

if "filter_otable_nosynteny" not in config:
    config["filter_otable_nosynteny"] = 'n'

### WORKFLOW

#Final output
rule Target:
    input: "SCORPiOs_"+config['jobname']+"/.cleanup_"+str(ITER)

#include the 5 modules
include: "module_build_trees.snake"
include: "module_orthology_table.snake"
include: "module_synteny_ortho_para.snake"
include: "module_graphs_orthogroups.snake"
include: "module_correct_trees.snake"
