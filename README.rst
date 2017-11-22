=======================
Coupled Moves Practical
=======================

Introduction
------------

In this activity, you will utilize the Coupled Moves [NO2015]_ protocol within Rosetta to find mutations capable of redesigning enzyme specificity for an alternate target ligand

Coupled Moves utilizes a flexible backbone modeling approach, achieving a better performance than a comparable fixed backbone method for 16 out of 17 tested benchmark cases [NO2015]_.

Coupled Moves operates according as follows (Fig 3. from [NO2015]_):

.. image:: journal.pcbi.1004335.g003.PNG
   :align: center
   :width: 70 %

Run Coupled Moves
-----------------

1. From within the coupled_moves .zip folder, open ``run.py`` in your text editor of choice.
#. Find the ``coupled_moves_path`` at the top of ``run.py`` and set it to the appropriate location of your compiled Rosetta coupled_moves binary.
#. Run ``python run.py``. The full command line call to each instance of Rosetta will be displayed, and will look something like this:

   ``/home/user/rosetta/rosetta_src_2017.45.59812_bundle/main/source/bin/coupled_moves.static.linuxgccrelease -s /home/coupled_moves/2O7B_HC4/2O7B_with_HC4.pdb -resfile /home/coupled_moves/2O7B_HC4/2O7B.resfile -extra_res_fa /home/coupled_moves/2O7B_HC4/HC4_from_2O7B.params -mute protocols.backrub.BackrubMover -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -coupled_moves::mc_kt 0.6 -coupled_moves::initial_repack false -coupled_moves::ligand_mode true -coupled_moves::fix_backbone false -coupled_moves::bias_sampling true -coupled_moves::boltzmann_kt 0.6 -coupled_moves::bump_check true -extra_res_fa /home/kyleb/algosb/coupled_moves/2O7B_HC4/MDO_from_2O7B.params``

   Important flags explained:

   * ``-ex1 -ex2 -extrachi_cutoff`` tell Rosetta's side chain packing algorithm to sample extra subrotamers for chi1 and chi2 angles of all side chains (`Packer documentation <https://www.rosettacommons.org/docs/latest/rosetta_basics/options/packing-options>`_).
   * ``-mute`` suppresses extraneous output from printing at the command line.
   * TODO

#. Output will be saved in a new directory named ``output``

Analysis
--------

Python package requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have provided full example output for the analysis stage. Run the analysis script as follows:

::

  python analyze_coupled_moves.py example_output/3KZO*

The analysis script will compare the distributions of output sequences all-by-all for all input output folders. Look for the lines that contain ``3KZO_SN0 over 3KZO_AN0``, which are mutations enriched in the non-native substrate (SN0/N-succinyl-L-ornithine) over the native substrate (AN0/N-acetyl-L-ornithine).

Additionally, a weblogo will be created in each output folder, and will look something like this:

|out1|  |out2|

.. |out1| image:: 3KZO_AN0-logo.pdf
   :width: 49 %
.. |out2| image:: 3KZO_SN0-logo.pdf
   :width: 49 %

Structure activity: load the known mutant crystal structure and compare one of the output structures.

.. [NO2015] Noah Ollikainen, René M. de Jong, and Tanja Kortemme. Coupling Protein Side-Chain and Backbone
   Flexibility Improves the Re-design of Protein-Ligand Specificity. *PLOS Comput Biol*, 11(9):e1004335,
   September 2015. ISSN 1553-7358. doi: 10.1371/journal.pcbi.1004335.
   URL http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004335

.. [CS2008] Colin A. Smith and Tanja Kortemme. Backrub-Like Backbone Simulation Recapitulates Natural Protein
   Conformational Variability and Improves Mutant Side-Chain Prediction. *Journal of Molecular Biology*, 380(4):
   742–756, July 2008. ISSN 0022-2836. doi: 10.1016/j.jmb.2008.05.023. URL http://www.sciencedirect.com/science/article/pii/S0022283608005779.
