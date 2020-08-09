# dmiGenerator

## Purpose
dmiGenerator generates a set of maximum theoretical* Dobzhansky and Muller incompatibilities given a genotype for N genes (e.g., ”AbCD”), if the user does not provided the genotype a random genotype for N genes where each gene is sampled as a random binomial to be either ancestral or derived. 

* All two-locus gene interactions that could occur in a perfect snowball scenario where one lineage evolves the given genotype and the other lineage evolves the exact opposite.

Examples:
```bash
dmiGenerator.py AbcDefG
dmiGenerator.py -n 16
dmiGenerator.py -v -n 26
dmiGenerator.py -v -g -n 6
dmiGenerator.py -vgn 6
dmiGenerator.py -t 18
```

## Options
`-h` shows program's help

`-v` activates verbose mode
```console
~$ dmiGenerator.py -v -n 6
Number of genes: 6
Original input (final genotype for lineage 1): ADEBCF
Final genotype for lineage 2): adebcf
Mutation per time in lineage 1: ['Adebcf', 'ADebcf', 'ADEbcf', 'ADEBcf', 'ADEBCf', 'ADEBCF']
Mutation per time in lineage 2: ['adebcf', 'adebcf', 'adebcf', 'adebcf', 'adebcf', 'adebcf']
Number of DMIs pairs: 15 (DI: 0, MI: 15)
Orr`s DMIs expection (K(K-1)/2): 15

Pair of genes with potential DMI:
['Da', 'Ea', 'Ed', 'Ba', 'Bd', 'Be', 'Ca', 'Cd', 'Ce', 'Cb', 'Fa', 'Fd', 'Fe', 'Fb', 'Fc']
```

`-t` replaces regular output for counts (# genes, # dmis, # ddis, # adis, # orrEstimate)
```console
~$ dmiGenerator.py -t -n 6
6       15      8       7       15.0
```

`-s` modifies output separator 
```console
~$ dmiGenerator.py -s ";" -n 6
Ab;Db;Da;Eb;Ea;Ed;FB;FA;FD;FE;CB;CA;CD;CE;Cf
```

`-g` shows a schematic representation of every mutation and its incompatibilities in both genotypes 
```console
~$ dmiGenerator.py -g -n 6
['Bd', 'Ed', 'Eb', 'CD', 'CB', 'CE', 'Fd', 'Fb', 'Fe', 'FC', 'AD', 'AB', 'AE', 'Ac', 'AF']
|d    D| |d    D| |d    D| |d    D| |d    D| |d    D|
|b    b| |b    B| |b    B| |b    B| |b    B| |b    B|
|e    e| |e    e| |e    E| |e    E| |e    E| |e    E| 
|c    c| |c    c| |c    c| |C    c| |C    c| |C    c|
|f    f| |f    f| |f    f| |f    f| |f    F| |f    F|
|a    a| |a    a| |a    a| |a    a| |a    a| |A    a|
```
