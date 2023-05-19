
Transcription Binding Predictor
===============================

A Python project that predicts whether a sequence of DNA will be bound by a specific transcription factor in a given condition. 

Given a set of sequences which contains a mix of sequences that are bound by the transcription factor or those that are unbound. The program predicts a subset of the sequences that are likely bound.

Given a PWM (positional weight matrix) previously discovered on known bound sequences to the transcription factor using ChipSeq, the program outputs the top 2000 sequences with the highest probability of being bound, taking account of reverse sequence complements as well. 

The output of the top 2000 most likely bound sequences will be output in a "predictions.txt" that looks like: 

seq14307
seq1929
seq10381
seq20681
seq10030
seq16200
seq8786
seq2876
seq13921
...

Sequences are ranked from the highest proability of binding to the lowest out of the 2000 top sequences.


Deliverables:
-------------

transcription_binding_predictor.py -- code for predicting bound sequences 

predictions.txt -- output of the top 2000 most likely bound sequences 

predictions.zip -- zipped csv of predictions.txt


Usage
-----

To run the program, navigate to the project directory and run:

> python3 transcription_binding_predictor.py test_reads.fasta pwm.txt 

The program takes the following arguments:

* `--test_reads.fasta`: A fasta file of the test reads that contain bound and unbound sequences.
* `--pwm.txt`: a PWM previously discovered on known sequences.

Examples
--------

Here is an example of how to run the program: 

>python3 transcription_binding_predictor.py test.fasta project2a_PWM.txt 


Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For ~25,000 test reads and a PWM of 21 by 4 the run time is:

real	1m22.459s
user	1m22.025s
sys	0m0.636s
