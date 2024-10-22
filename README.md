# A Quantitative Microbial Risk Assessment (QMRA) framework for exposure from toilet flushing using experimental aerosol concentration measurements

Ciara A. Higham<sup>1,*</sup>, Dr Martín López-García<sup>2</sup>, Prof. Catherine J. Noakes<sup>3</sup>, Emma Tidswell<sup>3</sup> and Dr Louise Fletcher<sup>3</sup>

<sup>1</sup> EPSRC Centre for Doctoral Training in Fluid Dynamics, University of Leeds, Woodhouse Lane, Leeds, LS2 9JT, United Kingdom.  
<sup>2</sup> School of Mathematics, University of Leeds, Woodhouse Lane, Leeds, LS2 9JT, United Kingdom.  
<sup>3</sup> School of Civil Engineering, University of Leeds, Woodhouse Lane, Leeds, LS2 9JT, United Kingdom.  
<sup>*</sup> Corresponding email: sccah@leeds.ac.uk

Supporting code and data for "A Quantitative Microbial Risk Assessment (QMRA) framework for exposure from toilet flushing using experimental aerosol concentration measurements".

The code is written in R version 4.3.2 and has been run using RStudio 2023.12.0.

## Repository contents
- `raw_data/`: a directory containing the raw particle concentration data as `.csv` files for
  - two scenarios: S1 vs S2
    - at particle counter locations: A and B
      - at 3 ventilation rates: 1.5, 3, 6 air changes per hour (ACH)
        - each with 3 replicates (a, b and c)
          
- `particle_conc/`: a directory containg the code to plot particle concentration as a time series
  -`particle_conc_plot.R`: R script that reads in the raw particle concentration data and plots a time series of the mean concentration across the three replicates. The cubicle scenario is plotted as a purple line and the no cubicle scenario is plotted as a green line. Standard errors are plotted as a shaded region. The particle sizes are split into further sub plots.
  
- `male_female_violin/`: a directory containing the code to model infection risk for $t_{male}$ and $t_{female}$ occupancy times
  - `dr_mc_male_female.R`: R script that reads in the raw particle concentration data and uses a stochastic Monte Carlo approach to quantify exposure and normalised infection risk to SARS-CoV-2 and norovirus for $t_{male}$ and $t_{female}$ occupancy time distributions, as defined in Table 1 and Fig. 2 of the supplementary material.
  - `infection_risk_violin.R`: R script that reads in the output data from running `dr_mc_male_female.R`, located in the `dr_output/` directory and plots the risk of infection (%) for SARS-CoV-2 and norovirus at particle counter location A and B for times $t = 0, 60, 240$ s.
  - `summary_stats.R`: R script that reads in the output data from running `dr_mc_male_female.R`, located in the `dr_output/` directory and gives the mean and median infection risks, details in Supplementary Table I.

- `heat_map/`: a directory containing the code to model infection risk for a uniform range of $t_{dur}$. Plots a heat map of normalised infection risk with varying $t_{enter}$ and $t_{dur}$.
  - `dr_mc_uniform_tdur.R`: R script that reads in the raw particle concentration data and uses a stochastic Monte Carlo approach to quantify exposure and infection risk to SARS-CoV-2 and norovirus for $t_{enter}$ [0-600] s and $t_{dur}$ [1-900] s in interval of 1 second.
  - `risk_heat_map.R`: R script that reads in the output data from running `dr_mc_uniform_tdur.R`, located in the `dr_output/` directory and plots heat map of normalised infection risk with varying $t_{enter}$ and $t_{dur}$.
  - `correlation.R`: R script that calculates the values of the Spearman correlation coefficients between variables and infection risk.

*Note: `particle_conc_plot.R` can be run independent of other scripts but `dr_mc_uniform_tdur.R` must be run prior to running `risk_heat_map.R` and `correlation.R`. `dr_mc_male_female.R` must be run prior to running `infection_risk_violin.R` and `summary_stats.R`.*

### Paper figures

- Fig. 4: `1.5_counter_a.png`, `3_counter_a.png`, `6_counter_a.png` located in`particle_conc/plots/counter_a/` after running `particle_conc_plot.R`.
- Fig. 5: `1.5_counter_b.png`, `3_counter_b.png`, `6_counter_b.png` located in`particle_conc/plots/counter_b/` after running `particle_conc_plot.R`.
- Fig. 6: `norovirus_counter_a.png` and  `sars_cov_2_counter_a.png` located in `male_female_violin/plots/counter_a/` after running `dr_mc_male_female.R` followed by `infection_risk_violin.R`.
- Fig. 7: `norovirus_S1_counter_a.png`, `norovirus_S2_counter_a.png` located in `heat_map/plots/norovirus` after running `dr_mc_uniform_tdur.R` followed by `risk_heat_map.R`.
- Fig. 8: `sars_cov_2_S1_counter_a.png`, `sars_cov_2_S2_counter_a.png` located in `heat_map/plots/sars_cov_2` after running `dr_mc_uniform_tdur.R` followed by `risk_heat_map.R`.




