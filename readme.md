# XPlode: Explaining Repaired Data with CFDs

> Based on the conference paper "Explaining Repaired Data with CFDs", by J. Rammelaere and F. Geerts, to be published in the Proceedings of the VLDB Endowment (PVLDB 2018). Code by Joeri Rammelaere.

A more interactive version of this code is available on the CodeOcean platform: https://codeocean.com/capsule/3719564/

## Quick Start Guide
The algorithm takes five arguments as input:
1. A dirty dataset
2. A (partially) cleaned version of the same dataset
3. Confidence leeway (referred to as Epsilon in the paper), between 0 and 1. CFDs will be discovered with a confidence of 1 minus this parameter
4. Support threshold (referred to as Delta in the paper), between 1 and the size of the dataset

## Available Datasets
The custex.csv and custdirty.csv files refer to the Customers dataset used as running example in the paper.
The other 4 dataset in the datasets/ folder are the clean versions of the datasets used in the experimental section.
In the dirtysets/ subfolder, for each dataset there are 3 dirtied datasets, with 10% errors generated for a single CFD.
The used CFDs can be found in the ${dataset}_cfds_order_x_.txt files. For the different experiments, the first _x_ cfds from these files were used.

Statistics of the datasets:

Dataset | Nr. Tuples | Nr. Attributes
:--- | :--- | :---
Abalone | 8354 | 9
Adult | 97684 | 11
Customers | 8 | 8
Soccer | 200000 | 10
SP500 | 245148 | 7

## Error Generation
We made use of the BART tool (http://www.db.unibas.it/projects/bart/) for introducing violations in the datasets. BART takes a
dataset and a set of data quality rules as input, and inserts
a predefined percentage of violations into the data. To get the required quality rules, we used our
implementation of CTane to discover 100% confident
CFDs on the datasets. Combining the clean
and dirty datasets, we generate partial repairs by starting
from the dirty dataset, and replacing a subset of the dirty
tuples with their clean variants. For each dataset, we obtain
3 different dirty datasets by using 3 different CFDs discovered
on the clean data, denoted by CFD 1, CFD 2, and
CFD 3. When considering a dirty dataset corresponding to
a CFD i, then CFD i is the target CFD. In other words, it is
the CFD we want to discover by repairing the corresponding dirty dataset. 
Obviously, CFD i is different for each dataset.

## Paper Abstract
Many popular data cleaning approaches are rule-based: 
Constraints are formulated in a logical framework, and data is considered dirty if constraints are violated.
These constraints are often discovered from data, but to ascertain their validity, user verification is necessary. 
Since the full set of discovered constraints is typically too large for manual inspection, recent research integrates user feedback into the discovery process. 
We propose a different approach that employs user interaction only at the start of the algorithm: 
a user manually cleans a small set of dirty tuples, and we infer the constraint underlying those repairs, called an explanation. 
We make use of conditional functional dependencies (CFDs) as the constraint formalism. 
We introduce XPlode, an on-demand algorithm which discovers the best explanation for a given repair. 
Guided by this explanation, data can then be cleaned using state-of-the-art CFD-based cleaning algorithms. 
Experiments on synthetic and real-world datasets show that the best explanation can typically be inferred using a limited number of modifications. 
Moreover, XPlode is substantially faster than discovering all CFDs that hold on a dataset, and is robust to noise in the modifications.
