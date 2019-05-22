# A Practical Guide To Estimating the Heritability of Pathogen Traits - Scripts and Data

Pathogen traits, such as the virulence of an infection, can vary significantly between patients. A major challenge is to measure the extent to which genetic differences between infecting strains explain the observed variation of the trait. This is quantified by the trait's broad-sense heritability, H<sup>2</sup>. A recent discrepancy between estimates of the heritability of HIV-virulence has opened a debate on the estimators' accuracy. Here, we show that the discrepancy originates from model limitations and important lifecycle differences between sexually reproducing organisms and transmittable pathogens. In particular, current quantitative genetics methods, such as donor-recipient regression (DR) of surveyed serodiscordant couples and the phylogenetic mixed model (PMM), are prone to underestimate H<sup>2</sup>, because they neglect or do not fit to the loss of resemblance between transmission partners caused by within-host evolution. In a phylogenetic analysis of 8,483 HIV patients from the UK, we show that the phenotypic correlation between transmission partners decays with the amount of within-host evolution of the virus. We reproduce this pattern in toy-model simulations and show that a phylogenetic Ornstein-Uhlenbeck model (POUMM) outperforms the PMM in capturing this correlation pattern and in quantifying H<sup>2</sup>. In particular, we show that POUMM outperforms PMM even in simulations without selection - as it captures the mentioned correlation pattern - which has not been appreciated until now. By cross-validating the POUMM estimates with ANOVA on closest phylogenetic pairs (ANOVA-CPP), we obtain H<sup>2</sup>~0.2, meaning about 20% of the variation in HIV-virulence is explained by the virus genome both for European and African data.

# License
The content of this project itself is licensed under the [Creative Commons Attribution 3.0](http://creativecommons.org/licenses/by/3.0/us/deed.en_US) license. The accompanying R-packages [toyepidemic](https://venelin.github.io/toyepidemic/index.html) and [patherit](https://venelin.github.io/patherit/index.html) are licensed under [version 3 of the GNU General Public License](http://www.gnu.org/licenses/).

# What is in this repository?

This repository contains the R-packages (not including 3rd party dependencies) as well as R-scripts used to perform toy-model simulations, analyze simulation results and analyze HIV-1 data for the article

*Venelin Mitov, Tanja Stadler; A Practical Guide to Estimating the Heritability of Pathogen Traits, Molecular Biology and Evolution, Volume 35, Issue 3, 1 March 2018, Pages 756–772, https://doi.org/10.1093/molbev/msx328*

All of this R-code was executed on R version 3.2.4 (2016-03-10) running on
x86_64-apple-darwin13.4.0 (64-bit). In addition, the toy-model simulations and ML fits of the POUMM and ANOVA-CPP methods were executed on the Brutus high-performance cluster at ETH Zurich (https://www1.ethz.ch/id/services/list/comp_zentral/cluster/index_EN).

## 1. Third party dependencies needed for running the scripts
In order to run the R-scripts you need to install the following dependencies:
xtable, knitr, rmarkdown, ape, coda, data.table, adaptMCMC and Rmpfr.

At the time of writing this documentation, all of the above packages were available on CRAN.

## 2. Installing in-house developed packages
You have to install the patherit and benchtable packages from their source-bundles included in
this directory. To do this, execute the following commands within R:

install.packages('patherit_0.4.tar.gz', repos=NULL, type='source')
install.packages('toyepidemic-0.1.0.tar.gz', repos=NULL, type='source')
install.packages('benchtable_1.0.tar.gz', repos=NULL, type='source')

NOTE: Included in the patherit package is the version of the POUMM package used by all R-scripts.
Since, the interface of the current version of the POUMM package has numerous changes w.r.t.
patherit v0.4, the R-scripts would fail to run with the newest POUMM version.

## 3. Simulations of the toy-model
To simulate the toy-model from scratch, you need to run individual chunks manually in the
SimulateToyModel.Rmd - file. This is a computationally expensinve and long-running task, which
needs the generation of benchtable shell-scripts, upload of these scripts together with
initialized data-files to a parallel cluster, execute parallel jobs, and gather result-files back
on a local computer, in order to analyze the results. You will need to install the packages patherit and
benchtable on the remote cluster as well. Read the instructions and R-code in the file
SimulateToyModel.Rmd (chunks initiate-data, simulate-and-fit-POUMM, collect-results).

## 4. Generating figures for the toy-model simulations

To generate the figures you first need to extract the archive-file 2018-MitovAndStadler-MBE-SimulationData.zip in the root directory of a working copy of this git repository. Due to its large file size (more than 7Gb), this file is not stored on github, but you can contact the corresponding author (Venelin Mitov) or the the last autor (Tanja Stadler) to get access to this file. This file should extracted in this directory prior to executing the chunks in this script. Then, all that needs to be done is executing the corresponding chunks in the file GenerateFigures.Rmd.

## 5. Analysis of HIV-data

To analyse the HIV-data, currently, you need a permission for accessing the data
from the UK drug resistance database. I recommend contacting the database at
mrcctu.hivrdb@ucl.ac.uk (CC: venelin.mitov@bsse.ethz.ch) and request the phylogeny
and spVL data published in Hodcroft et al. 2014 (see article main text for a full
reference). Once you obtain the permission from the UK-drug resistance database, 
please, write an e-mail to vmitov@gmail.com or tanja.stadler@bsse.ethz.ch. Then, 
I or Tanja will send you a a file called HIV.zip that should be extracted in the 
subdirectory DATA of this directory.

Then, you will be able to execute the commands in the R-script AnalyseHIVspVL.Rmd, 
as well as generate the HIV-data figures in GenerateFigures.Rmd.

For questions, please write to vmitov@gmail.com or to tanja.stadler@bsse.ethz.ch.
