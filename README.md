# Experiment Materials
This repository stores the experiments materials from the paper "A Meta-Learning
Algorithm Selection Approach for the Quadratic Assignemnt Problem" to be
pusblished at the Proceedings of [IEEE Congress on Evolutionary Computation
(IEEE CEC) 2018](http://www.ecomp.poli.br/~wcci2018/).

The content of each folder is described next.

### Optimization Algorithms
This folder contains the source-codes of the meta-heuristics with their default
parameters. It is also provided the set of seeds used for the experiments and
the runner scripts for the landscape features extraction.

### Datasets
These are the generated classification datasets. The full and cleansed versions
are given, along with the datasets used in the final cascade scheme.

### Cascade Reults
It contains the results metrics, such as the Accuracy and F-Scores, of the
cascaded model given by the script ```cascade_classifier.py```, which performs a
KFold Cross Validation using the best sets of features for each cascade level.

### QAP Intances
These are all the instances retrieved from
[QAPLIB](http://anjos.mgi.polymtl.ca/qaplib/inst.html).
