#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

"""
Generate from a given secuence of genes (or a random sequence) all potential DMIs

This version generates thousands of genes

The idea is produce a list of dmis in the following format:

0/1:2/1-1/1:2/0-0/1:2/1 == AB-Bc-AC

"""

import argparse
import numpy as np

#get some arguments
parser = argparse.ArgumentParser(description="Script find the potential DMIs given a sequence of loci")
parser.add_argument("-n", "--number",  metavar='N', type=int, default=3, dest="number", help="Generate a random sequence of N genes")
parser.add_argument("-s", "--seed", dest="seed", metavar='N', type=int, help="Seed")
parser.add_argument("-l", "--letters", dest="letters", help="Return in DMIs as letters (A-Z)", action="store_true")

#put them into args
args = parser.parse_args()

#create rng generator
rng = np.random.default_rng(args.seed)

# for thousands only use random generator
lineage1_new = []

#loop for every asked gene
for i in range(args.number):
    #define if current gene is ancestral (true) or derived (false)
    if rng.choice([True, False]):
        lineage1_new.append(f"{i}/0")
    else:
        lineage1_new.append(f"{i}/1")

# reorder each gene, to avoid 1 is always first 
rng.shuffle(lineage1_new)


### Some functions 

def _get_gene_state(gene_state):
    """Function for controling splits in one place"""
    return gene_state.split("/")


def _isderived(gene):
    """Replace of isupper check in thousands versions
    This verify if the given gene is derived or not.
    """    
    _gene, _state = gene.split("/")
    return _state == "1"


def _convert_numbers_to_letters(gene):
    """Convert gene in number format to letters (A-Z) being lower case ancestral
    """
    gene, state = _get_gene_state(gene)
    if int(gene) <= 25:
        letter = chr(65 + int(gene))
        if state == "0":
            letter = letter.lower() 
        return letter


def reconstruct_lineage_history(lineage):
    """Reconstruct the history of mutation (in a finite site model) of every gene
        A genome history ABC will have three slides in the history (ordered from past to present): Abc, ABc, ABC
    """
    lineage_new_history = []
    for index in range(len(lineage)):
        slice = []
        #split the genes increasing one by one, to mantain this part in the
        #same state as the lineage (present time), in the first round it a lineage
        #with this genome history:  AbCDEf will be Abcdef
        for _genestate in lineage[:index + 1]:
            slice.append(_genestate)

        # the rest of genes (in this case decreasing one by one), are change to
        # the anscestral state
        for _genestate in lineage[index + 1:]:
            _gene, _state = _get_gene_state(_genestate)
            slice.append(f"{_gene}/0")
        
        # append every slice 
        lineage_new_history.append(slice)

    return lineage_new_history

# create the hypotetical lineage 2 as a oposite mirror of lineage 1
lineage2_new = []
for i in lineage1_new:
    _gene, _state =  _get_gene_state(i)
    _opposite_state = 1 - int(_state)
    lineage2_new.append(f"{_gene}/{_opposite_state}")



# get histories for each lineagef
lineage1_new_history = reconstruct_lineage_history(lineage1_new)
lineage2_new_history = reconstruct_lineage_history(lineage2_new)


dmis_new = []


# print(f"{lineage1_new_history=}")
# print(f"{lineage2_new_history=}")



##thousands version
#do a loop in the first lineage history, saving the index.
#do a loop in the first genotype, in the first letter, check if have or not a mutation
    #if has the mutation iterate with all genes in lineage 2 in that time (or index)
    #if does not have a mutation move to the lineage 2 and does the iteration with all genes in lineage 1
#while is getting all possible combination, do the evauation, if Aa noDMI, if pair is in the genotype is not a DMI.
#But when the pair is not in the genotype is a potential DMI
for index, genotype_l1_new in enumerate(lineage1_new_history):
    # print(f"{genotype_l1_new=}")
    # print(f"{genotype_l1_new[index]=}")
    gene_l1_new = genotype_l1_new[index]
    if _isderived(gene_l1_new):
        # print(f"{genotype_l1_new[index]=}")
        for comparison_new in lineage2_new_history[index]:
            # print(f"{comparison_new=}")
            pair_new = [gene_l1_new, comparison_new]
            # print(f"{pair_new=}")
            # print([c in pair_new for c in genotype_l1_new], pair_new, genotype_l1_new)
            if sum([c in pair_new for c in genotype_l1_new]) < len(pair_new) and _get_gene_state(gene_l1_new)[0] != _get_gene_state(comparison_new)[0]:
                dmis_new.append(pair_new)
    else:
        gene_l2_new = lineage2_new_history[index][index]
        # print(f"{gene_l2_new=}")
        if _isderived(gene_l2_new):
            for comparison_new in genotype_l1_new:
                pair_new = [gene_l2_new, comparison_new]
                if sum([c in pair_new for c in lineage2_new_history[index]]) < len(pair_new) and _get_gene_state(gene_l2_new)[0] != _get_gene_state(comparison_new)[0]:
                    dmis_new.append(pair_new)


# output
#convert number notation to letters
formated_output = []
if args.letters:
    if args.number <= 26:        
        for d in dmis_new:
            formated_dmi = ""
            for g in d:
                formated_dmi += _convert_numbers_to_letters(g)     
            formated_output.append(formated_dmi)   
        print("-".join(formated_output))
    else:
        print("There are only 26 letters in the roman alphabet. Reduce the number of genes to maximum 26")
else:
    for d in dmis_new:
        formated_output.append(":".join(d))   
    print("-".join(formated_output))



