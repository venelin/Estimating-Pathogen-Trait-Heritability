# A Practical Guide To Estimating the Heritability of Pathogen Traits - Scripts for the Simulations and the Analysis of Real Data

Pathogen traits, such as the virulence of an infection, can vary significantly between patients. A major challenge is to measure the extent to which genetic differences between infecting strains explain the observed variation of the trait. This is quantified by the trait's broad-sense heritability, $H^2$. A recent discrepancy between estimates of the heritability of HIV-virulence has opened a debate on the estimators' accuracy. Here, we show that the discrepancy originates from model limitations and important lifecycle differences between sexually reproducing organisms and transmittable pathogens. In particular, current quantitative genetics methods, such as donor-recipient regression (DR) of surveyed serodiscordant couples and the phylogenetic mixed model (PMM), are prone to underestimate $H^2$, because they neglect or do not fit to the loss of resemblance between transmission partners caused by within-host evolution. In a phylogenetic analysis of 8,483 HIV patients from the UK, we show that the phenotypic correlation between transmission partners decays with the amount of within-host evolution of the virus. We reproduce this pattern in toy-model simulations and show that a phylogenetic Ornstein-Uhlenbeck model (POUMM) outperforms the PMM in capturing this correlation pattern and in quantifying $H^2$. In particular, we show that POUMM outperforms PMM even in simulations without selection - as it captures the mentioned correlation pattern - which has not been appreciated until now. By cross-validating the POUMM estimates with ANOVA on closest phylogenetic pairs (ANOVA-CPP), we obtain $H^2\approx 0.2$, meaning about 20\% of the variation in HIV-virulence is explained by the virus genome both for European and African data.

# License
The content of this project itself is licensed under the [Creative Commons Attribution 3.0](http://creativecommons.org/licenses/by/3.0/us/deed.en_US) license. The accompanying R-packages are licensed under [version 3 of the GNU General Public License](http://www.gnu.org/licenses/).
