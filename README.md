# maximum_power_info_ratchet
Scripts and data for manuscript "Maximizing power and velocity of an information engine" by Tushar K. Saha, Joseph N. E. Lucero, Jannik Ehrich, David A. Sivak, and John Bechhoefer

## Directories reference

data/: Directory containing data to generate all figures in main text and SI appendix.

doc/: Directory containing figures found in the main text and SI appendix. 

src/example_exp_analysis/: Directory containing a self-contained example of how analysis of experimental data is done. README file contained in the directory explains further. 

src/analysis.ipynb: Jupyter notebook containing theoretical calculations that regenerate the curves found in the main text. Also has details about the calculations that can be found in the SI Appendix. 

src/propagators.nb: Mathematica notebook containing the details of how the propagators, discussed in the SI Appendix section E and N, are derived. Expressions derived are used in the Jupyter notebook above.

src/quadratic-to-linear_potential/: Directory containing script carrying out the calculations detailed in SI Appendix Section K.1 and generates Figs. S7A and S7B. Script outputs data file to data/ directory and figure to doc/ directory.

src/noisy_analysis/: Directory containing scripts (and corresponding derivation in the Mathematica notebook) to compute the efficiency of the engine given empirical measurement noise, as detailed in SI Appendix Section N.
