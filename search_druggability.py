import sys
import re
from collections import defaultdict


#Input: Drug List
#Output: Drugs, synonyms, side effects and indications

#For each drug in the drug list, the script finds all synonym names. For each drug name + synonyms, the script searches for indications and side effects.

#USAGE:
#python search_druggability.py chemical.aliases.v3.1.tsv adverse_effects_raw.tsv indications_raw.tsv label_mapping.tsv druglist.txt output.txt

try:
    aliases = open(sys.argv[1], "r")
    adverse = open(sys.argv[2], "r")
    indications = open(sys.argv[3], "r")
    label_mapping = open(sys.argv[4], "r")
    drugs = open(sys.argv[5], "r")
    output = open(sys.argv[6], "w")
except:
    sys.exit("File does not exist!")


#CREATES DICTIONARY OF SYNONYMS (key=cid; value=synonyms)
synonyms = defaultdict(list) #creates dictionary of chemical aliases
for line in aliases: #reads file of the names of medicines
    splitline=line.split("\t")
    cid = splitline[0].upper() #assigns cid to values of the first column of the file
    name = splitline[1].upper() #assigns name to values of the second column of the file
    synonyms[cid].append(name) #assigns newcid2 as keys and names as values
aliases.close()

#CREATES DICTIONARY OF LABEL_MAPPING (key=drug; value=label_ids)
label_mapping_dict = defaultdict(list) #creates dictionary of label mapping
for line in label_mapping: #reads file of label mapping
    splitline=line.split("\t")
    label_name = splitline[1].upper() #assigns names to values of the second column of the file
    label_mapping_id = splitline[6].rstrip().upper() #assigns label_mapping_id to values of the seventh column of the file
    label_mapping_dict[label_name].append(label_mapping_id) #assigns values to the dictionary
label_mapping.close()

#CREATES DICTIONARY OF ADVERSE EFFECTS (key=label_id; values=adverse effects)
adverse_dict = defaultdict(list) #creates dictionary of adverse effects
for line in adverse: #reads file of adverse effects
    splitline=line.split("\t")
    label_id = splitline[0].upper() #assigns label_id to values of the first column of the file
    effect = splitline[2].rstrip().upper() #assigns side effect to values of the third column of the file
    adverse_dict[label_id].append(effect)
adverse.close()

#CREATES A DICTIONARY OF ALL ADVERSE EFFECTS (key=drug;values=all adverse effects)
colaterals = defaultdict(list)
for key, value in label_mapping_dict.iteritems():
    drugname = key
    for code in label_mapping_dict.get(key):
        for key,value in adverse_dict.iteritems():
            if code == key:
                for i in adverse_dict.get(key):
                    if i not in colaterals[drugname]:
                        colaterals[drugname].append(i)
adverse_dict.clear()                      

#CREATES DICTIONARY OF INDICATIONS (key=label_id; values=indications)
ind = defaultdict(list)
for line in indications: #reads file of indications
    splitline=line.split("\t")
    label_id2 = splitline[0].upper() #assigns label_id2 to values of the first column of the file
    indication = splitline[2].rstrip().upper() #assigns indications to values of the third column of the file
    ind[label_id2].append(indication) #appends values of side effects to the respective key
indications.close()

#CREATES A DICTIONARY OF ALL INDICATIONS (key=drug; values= all indications)
indications_all = defaultdict(list)
for key, value in label_mapping_dict.iteritems():
    drugname = key
    for code in label_mapping_dict.get(key):
        for key,value in ind.iteritems():
            if code == key:
                for i in ind.get(key):
                    if i not in indications_all[drugname]:
                        indications_all[drugname].append(i)
ind.clear()

druglist = drugs.read().splitlines()
for drug in druglist:
    drug = drug.upper()
    print drug
    for key, value in synonyms.iteritems():
        if drug in synonyms.get(key):
            drug_synonyms = synonyms.get(key)
            inds = set()
            adv_eff = set()
            written = []
            for drug_synonym in drug_synonyms:
                for key,value in indications_all.iteritems():
                    if drug_synonym == key:
                        for indication in indications_all.get(key):
                            inds.add(indication)
                for key, value in colaterals.iteritems():
                    if drug == key:
                        for efeito in colaterals.get(key):
                            adv_eff.add(efeito)
            if len(written)==0 or drug not in written:
                output.write(drug + '\t' + ','.join(drug_synonyms) + '\t' + ','.join(inds) + '\t' + ','.join(adv_eff) + "\n\n")
                print drug + '\t' + ','.join(drug_synonyms) + '\t' + ','.join(inds) + '\t' + ','.join(adv_eff) + "\n\n"
                written.append(drug)
output.close()
