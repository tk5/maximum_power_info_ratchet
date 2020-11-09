# maximum_power_info_ratchet
Scripts and data for manuscript "Maximizing power and velocity of an information engine" by Tushar K. Saha, Joseph N. E. Lucero, Jannik Ehrich, David A. Sivak, and John Bechhoefer

## Directories reference

data/: Directory containing data to generate all figures in main text and SI appendix.

doc/: Directory containing figures found in the main text and SI appendix. 

src/experimental_analysis/: Directory containing IGOR files that perform analysis on the experimental data and generate the associated figures in main text and SI appendix. 

src/experimental_analysis/example_analysis/: Directory containing a self-contained example of how analysis of experimental data is done. README file contained in the directory explains further. 

src/limits.py: Script computing expected limits of performance following discussion detailed in SI Appendix Section L. 

src/power_vs_frequency_numeric/: Directory containing script corresponding to the calculation detailed in SI Appendix Section E and generates the theoretical curve found in Fig. 3A of the main text. Script outputs to data/ directory.

src/power_vs_threshold_numeric/: Directory containing script corresponding to calculation that generates the theoretical curve found in Fig. 3B of the main text. Script outputs to data/ directory.

src/quadratic-to-linear_potential/: Directory containing script carrying out the calculations detailed in SI Appendix Section K.1 and generates Figs. S7A and S7B. Script outputs data file to data/ directory and figure to doc/ directory.
