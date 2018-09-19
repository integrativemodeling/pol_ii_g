   	
"""This script samples RNA_POL2_HUMAN subunits with GDOWN1 inhibitor cross-links data from CHAIT's GROUP
"""
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.mmcif
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.crosslinking
import os,sys

def add_atomic_rep(mol,chain,unstructured_bead_size,clr,prot):
        atomic = mol.add_structure('../data/'+prot +".pdb",chain_id=chain,offset=0)
        mol.add_representation(atomic, resolutions=[1,10],color = clr)
        if len(mol[:]-atomic) >0:
            mol.add_representation(mol[:]-atomic,resolutions=[10],color=clr)
        return mol

def add_nonatomic_rep(mol,chain,unstructured_bead_size,clr):
        mol.add_representation(mol[:],resolutions=[unstructured_bead_size],color=clr,setup_particles_as_densities=True,density_force_compute=False)
        return mol

##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
###################### SYSTEM SETUP #####################
num_frames = 50000
if '--test' in sys.argv:
    num_frames=25
    taskid=0
elif '--dry-run' in sys.argv:
    taskid=0
else:
    taskid = int(sys.argv[1])
#input files 
xlink_file = '../data/allcxg.csv'
all_cg_bead_size = 10

#MC parameters
RB_MAX_TRANS = 4.0
RB_MAX_ROT = 1.0
FLEX_MAX_TRANS = 4.0
SRB_MAX_TRANS = 1.0
SRB_MAX_ROT = 0.1


# Input sequences
pol2g_seq_file ='pol2gdown.fasta'
pol2g_seqs = IMP.pmi.topology.Sequences('../data/'+pol2g_seq_file)
pol2g_components={"A":["RPB1"],"B":["RPB2"],"C":["RPB3"],"D":["RPB4"],"E":["RPB5"],"F":["RPB6"],"G":["RPB7"],"H":["RPB8"],"I":["RPB9"],"J":["RPB10"],"K":["RPB11"],"L":["RPB12"],"X":["GDOWN1"]}
pol2g_colors ={"A":[0.1],"B":[0.2],"C":[0.3],"D":[0.4],"E":[0.5],"F":[0.6],"G":[0.7],"H":[0.8],"I":[0.9],"J":[0.05],"K":[0.15],"L":[0.25],"X":[0.86]}


# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)

if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('pol_ii_g.cif', 'w'))
    s.add_protocol_output(po)
    po.system.title = ('Architecture of Pol II(G) and molecular mechanism '
                       'of transcription regulation by Gdown1')

s.dry_run = '--dry-run' in sys.argv

st = s.create_state()

# Add Molecules for each component as well as representations
mols = []
pol2g_mols = []
gdown1=[]
pol=[]
#loop over for single copy of proteins with PDB structures
for chain in ['A','B','C','D','E','F','G','H','I','J','K','L']:
        prot = pol2g_components[chain][0]
        color = pol2g_colors[chain][0]
        mol1 = st.create_molecule(prot, sequence=pol2g_seqs[prot+chain],chain_id=chain)
        mol = add_atomic_rep(mol1,chain,all_cg_bead_size,color,'5flm')
       	mols.append(mol)
        pol2g_mols.append(mol)
        pol.append(mol)
for chain in ['X']:
        prot = pol2g_components[chain][0]
        color = pol2g_colors[chain][0]
        mol1 = st.create_molecule(prot, sequence=pol2g_seqs[prot+chain],chain_id=chain)
        mol = add_nonatomic_rep(mol1,chain,all_cg_bead_size,color)
        mols.append(mol)
        pol2g_mols.append(mol)
        gdown1.append(mol)
##  calling System.build() creates all States and Molecules (and their representations)
##  Once you call build(), anything without representation is destroyed.
##  You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
##  However these functions will only return BUILT representations
root_hier = s.build()
## Uncomment this for verbose output of the representation
IMP.atom.show_with_representations(root_hier)
## output to RMF
fname = 'test.rmf'
rh = RMF.create_rmf_file(fname)
IMP.rmf.add_hierarchy(rh, root_hier)
IMP.rmf.save_frame(rh)
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
# Setup degrees of freedom
#  The DOF functions automatically select all resolutions
#  Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

# DOF for POL2G
# Each structured unit is a single rigid body, with flexible beads corresponding to missing regions, and each subunit is a super rigid body
#swisnf_structured=[]
#swisnf_unstructured=[]
pol2g_complex=[]
pol2g_complexn=[]


for i,mol in enumerate(pol2g_mols):
  # create a rigid body for each domain with structural information, the flexible beads inside are part of the rigid body, the remaining part of the 
 #crystal structures of RPB1 and RPB2   
        if i in range(12):
		print i, mol
                pol2g_unstructured=mol.get_non_atomic_residues()
                pol2g_structured=mol.get_atomic_residues()
                pol2g_all=mol.get_residues()
                pol2g_complex.append(pol2g_all)
                pol2g_complexn.append(pol2g_unstructured)
     
	if i in [12]:
                pol2g_u=mol.get_non_atomic_residues()
                pol2g_rigid1=mol.get_residues()
                dof.create_flexible_beads(pol2g_rigid1,max_trans=FLEX_MAX_TRANS,resolution=10)


print dof

pol2g_molc2=[x for x in pol2g_complex]
pol2g_molc3=[x for x in pol2g_complexn]
dof.create_flexible_beads(pol2g_molc3,max_trans=FLEX_MAX_TRANS,resolution=10)


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################    
####################### RESTRAINTS #####################
output_objects = [] # keep a list of functions that need to be reported
display_restraints = [] # display as springs in RMF
 
# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
for mol in pol2g_mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    output_objects.append(cr)
    crs.append(cr)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print sf.evaluate(False)

# Excluded volume - automatically more efficient due to rigid bodies

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = gdown1,resolution=10)
ev1.add_to_model()
ev1.set_label('gdown1')
output_objects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print "ev1", sf.evaluate(False)

ev2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = gdown1,
                                                              other_objects = pol,
                                                              resolution=10)

ev2.add_to_model()
ev2.rs.set_weight(1.0)
ev2.set_label('all')
output_objects.append(ev2)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print "ev2", sf.evaluate(False)
print gdown1
print pol


ev3 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = pol2g_molc3,resolution=10)

ev3.add_to_model()
ev3.rs.set_weight(1.0)
ev3.set_label('all')
output_objects.append(ev3)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print "ev3", sf.evaluate(False)


print pol2g_mols
## Crosslink restraint
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb.create_set_from_file(xlink_file)

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,CrossLinkDataBase=xldb,length=21.0,label="XLDSS",filelabel='DSS',resolution=1,slope=0.01)
xlr.add_to_model()
output_objects.append(xlr)
display_restraints.append(xlr)
xlr.rs.set_weight(50.0)
# Set PSI value instead of sampling it
xlr.set_psi_is_sampled(True)
#psi = xlr.psi_dictionary["PSI"][0]
#psi.set_scale(0.04) 
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print sf.evaluate(False)

dof.get_nuisances_from_restraint(xlr) # needed to sample the nuisance particles (noise params)

####################### SAMPLING #####################
# First shuffle the system

if not s.dry_run:
    IMP.pmi.tools.shuffle_configuration(gdown1)

    dof.optimize_flexible_beads(100)
print dof.get_movers()
print len(dof.get_movers())
print gdown1
# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    crosslink_restraints=display_restraints,                     # will display like XLs
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
				    replica_exchange_maximum_temperature = 2.5,
				    num_sample_rounds = 1,
				    number_of_best_scoring_models = 5,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='gdownrb_%d' % taskid,
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_frames,
                                    test_mode=s.dry_run)

rex.execute_macro()


if '--mmcif' in sys.argv:
    po.flush()
