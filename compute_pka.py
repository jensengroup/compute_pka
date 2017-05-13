import sys
import math
import subprocess
sys.path.append("/usr/local/lib/python2.7/site-packages/")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import small_ref_set as rf

filename = sys.argv[1]
filename2 = sys.argv[2]
energies_file = open(filename, "r")
smiles_file = open(filename2, "r")

filename3 = filename.split(".")[0]
filename3 = "reference_"+filename3.split("_")[1]+"_"+filename3.split("_")[2]+".energies"

reference_energies_file = open(filename3, "r")

energies = {}
log_files = {}
names = []
ions = []
ion2smiles = {}
ref_energies = {}
atom_radius = 2
micro_pka_cutoff = 0.6

smarts_ref = ( ('[NX4+]',10.4,'N-sp3'),
               ('[NX3+]',10.9,'N-sp2'),
               ('[nX3+]',5.2,'N-heterocycle'),
               ('[OX2H1]',10.0,'O-sp2'),
               ('[SX2H1]',10.6,'S') )

def find_ref(reactant_name,product_name):

    reactant_smiles = ion2smiles[reactant_name] 
    product_smiles = ion2smiles[product_name] 

    reactant = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)

    r_smiles_input = "-:"+reactant_smiles
    reactant_inchi =  subprocess.check_output(['obabel','---errorlevel','0',r_smiles_input,'-oinchi','-xF']).split()[0]
    reactant = Chem.MolFromInchi(reactant_inchi)
    p_smiles_input = "-:"+product_smiles
    product_inchi =  subprocess.check_output(['obabel','---errorlevel','0',p_smiles_input,'-oinchi','-xF']).split()[0]
    product = Chem.MolFromInchi(product_inchi)
#    print Chem.MolToSmiles(reactant),Chem.MolToSmiles(product)
    
#    deprotonated_atom = fda.find_deprotonated_atom(reactant,product)

    charged_atoms = Chem.MolFromSmarts('[+,OX2H1,SX2H1]')
# find which atom is deprotonated
    charged_atoms_in_reactant = reactant.GetSubstructMatches(charged_atoms)
    charged_atoms_in_product = product.GetSubstructMatches(charged_atoms)

#    print "charged_atoms_in_reactant", charged_atoms_in_reactant,reactant_inchi
#    print "charged_atoms_in_product", charged_atoms_in_product, product_inchi

    deprotonated_atoms = [x[0] for x in charged_atoms_in_reactant if x not in charged_atoms_in_product]
    if len(deprotonated_atoms) == 0:
       print "hello"
       deprotonated_atoms = [x[0] for x in charged_atoms_in_product if x not in charged_atoms_in_reactant]

#    print "deprotonated_atoms",deprotonated_atoms
    if len(deprotonated_atoms) == 1: 
       deprotonated_atom = deprotonated_atoms[0]    
    elif len(deprotonated_atoms) > 1:
       deprotonated_atom = -1
    
# -1 means that more than one protonation state changes on ging from reactants to products

    if deprotonated_atom == -1:
        nref = "skip"
        pKa_nref = 0.0
#        print "warning", deprotonated_atom, reactant_smiles, product_smiles
        return nref, deprotonated_atom, pKa_nref, reactant_smiles, product_smiles;
  
    max_length = 0
    nref = "none"
    pKa_nref = 0.0
    for smarts, pKa, name in smarts_ref:
        ref_mol = Chem.MolFromSmarts(smarts)

        matches_tuple = reactant.GetSubstructMatches(ref_mol) 
        matches_list = [element for tupl in matches_tuple for element in tupl]
#        print smarts, deprotonated_atom,matches_list
        if reactant.HasSubstructMatch(ref_mol) and deprotonated_atom in matches_list:
#          print smarts, deprotonated_atom,list(reactant.GetSubstructMatches(ref_mol))
           length = len(list(reactant.GetSubstructMatches(ref_mol)))
           if length > max_length:
              nref = name
              max_length = length
              pKa_nref = pKa

#    if nref == "none": print "warning", reactant_name, product_name, reactant_smiles, product_smiles

    if nref == "N-sp3":
       nref = "ethylamine_1"
       better_ref, better_pKa_ref = find_better_ref(reactant, deprotonated_atom, rf.N_sp3_smiles)
       if better_ref != "none":
           nref = better_ref 
           pKa_nref = better_pKa_ref
       
    elif nref == "N-sp2":
       nref = "ethanimine_1" 
       better_ref, better_pKa_ref = find_better_ref(reactant, deprotonated_atom, rf.N_sp2_smiles)
       if better_ref != "none":
           nref = better_ref 
           pKa_nref = better_pKa_ref
              
    elif nref == "N-heterocycle":
       nref = "pyridine_1" 
       better_ref, better_pKa_ref = find_better_ref(reactant, deprotonated_atom, rf.N_heterocycle_smiles)
       if better_ref != "none":
           nref = better_ref 
           pKa_nref = better_pKa_ref
              
    elif nref == "O-sp2":
       nref = "phenol_0" 
       better_ref, better_pKa_ref = find_better_ref(reactant, deprotonated_atom, rf.O_sp2_smiles)
       if better_ref != "none":
           nref = better_ref 
           pKa_nref = better_pKa_ref

    elif nref == "S":
       nref = "ethanethiol_0" 
       better_ref, better_pKa_ref = find_better_ref(reactant, deprotonated_atom, rf.S_smiles)
       if better_ref != "none":
           nref = better_ref 
           pKa_nref = better_pKa_ref
              
    return nref, deprotonated_atom, pKa_nref, reactant_smiles, product_smiles;

def find_better_ref(reactant,deprot_atom, ref_list):
    nref = "none"
    pKa_nref = 0.0
    max_length = 0
    for smiles, pKa, name in ref_list:
        ref_mol = Chem.MolFromSmiles(smiles)
#       ref_mol = Chem.MolFromSmarts(smiles)
 
    #    print Chem.MolToSmiles(reactant), smiles,deprot_atom, list(reactant.GetSubstructMatch(ref_mol))

        if reactant.HasSubstructMatch(ref_mol) and deprot_atom in list(reactant.GetSubstructMatch(ref_mol)):
           length = len(list(reactant.GetSubstructMatch(ref_mol)))
           if length > max_length:
              nref = name
              max_length = length
              pKa_nref = pKa

    return nref, pKa_nref

for line in reference_energies_file:
    words = line.split()
    log_file  = words[0]
    ion = words[0].split("=")[0]
    if words[1] == "CURRENT":
        energy = float(words[-1])
    else:
        energy = float(words[6])
    if ion not in ref_energies:
       ref_energies[ion] = energy
    elif energy < ref_energies[ion]:
       ref_energies[ion] = energy

for line in smiles_file:
    words = line.split()
    ion = words[0]
    smiles = words[1]
    ion2smiles[ion] = smiles 

for line in energies_file:
    words = line.split()
    log_file  = words[0]
    ion = words[0].split("+")[0]
    name = ion.split("=")[0]
    if words[1] == "CURRENT":
        energy = float(words[-1])
    else:
        energy = float(words[6])
    if name not in names:
        names.append(name)
    if ion not in ions:
        ions.append(ion)
    if ion not in energies:
       energies[ion] = energy
       log_files[ion] = log_file
    elif energy < energies[ion]:
       energies[ion] = energy
       log_files[ion] = log_file

for name in names:
    charge = int(name.split("_")[1])    
    name_of_deprotonated = name.split("_")[0]+"_"+str(charge-1)
    if name_of_deprotonated not in names:
        continue

    protonated = [item for item in ions if name in item]
    deprotonated = [item for item in ions if name_of_deprotonated in item]

    sum_deprot = 0.0
    min_deprot_pka = 99999.9
    pka_dict = {}
    ref_dict = {}
    deprot_dict = {}
    for deprot in deprotonated:
        sum_prot = 0.0
        max_prot_pka = -99999.9
        prot_pka = []
        for prot in protonated:
            ref,deprotonated_atom, pKa_ref, r_smiles, p_smiles = find_ref(prot,deprot)
            if ref == "skip": 
                continue
            delta_G = energies[prot] - energies[deprot]
            ref_deprot = ref.split("_")[0]+"_"+str(int(ref.split("_")[1])-1)
            delta_G_ref = ref_energies[ref] - ref_energies[ref_deprot]
            pKa = pKa_ref + (delta_G_ref - delta_G)/1.36
            sum_prot += 10**pKa
#           if name == "Clozapine_2": print "Clozapine_2", prot,deprot,ref, pKa, pKa_ref, energies[prot] , energies[deprot],ref_energies[ref], ref_energies[ref_deprot],sum_prot
            index = prot+":"+deprot
            pka_dict[index] = pKa
            ref_dict[index] = ref.split("_")[0]+"_"+str(deprotonated_atom)
            prot_pka.append(pKa)
            if pKa > max_prot_pka: 
                max_prot_pka = pKa  


#       if ref != "skip": 
        if sum_prot > 0.0:
           sum_deprot += 1.0/sum_prot
           if max_prot_pka < min_deprot_pka: 
              min_deprot_pka = max_prot_pka
           deprot_dict[deprot] = prot_pka
#          if name == "Clozapine_2": print "Clozapine_2",prot,deprot,deprot_dict[deprot],sum_deprot

    ref_list = []
    ref_list_key = []
    for deprot in deprot_dict:
        micro_pka_list = []
        if abs(max(deprot_dict[deprot]) - min_deprot_pka) < micro_pka_cutoff:
            for micro_pka in deprot_dict[deprot]:
               if abs(micro_pka - min_deprot_pka) < micro_pka_cutoff:
                   micro_pka_list.append(micro_pka)

        for micro_pka in micro_pka_list:
            for key,pka_in_dict in pka_dict.iteritems():
               if pka_in_dict == micro_pka:
                   if ref_dict[key] not in ref_list:
                       ref_list.append(ref_dict[key])
                   if key not in ref_list_key:
                       ref_list_key.append(key)

#       if name == "Clozapine_2": print "Clozapine_2", micro_pka_list, ref_list, ref_list_key
    references = ""
    for reference in ref_list:
        references += ","+reference
    references = references[1:]

    pKa_app = -math.log(sum_deprot,10)

    low_energies = []
    for key in ref_list_key:
        prot = key.split(":")[0]
        deprot = key.split(":")[1]
        low_energies.append(energies[prot])
        low_energies.append(energies[deprot])

    print name, pKa_app, references, low_energies
