
#===============  P-value Simulation Procedures:  ========================

We conducted simulation to assess the effect of model misspecification on inference for the mediation edge (the edge between the two genes in a trio).

- We perform simulations based on real trios from WholeBlood tissue in GTEx

- We only consider trios that were significant under the GMAC inference as candidates for simulation (508 total trios for   Whole Blood). in the GMAC_postproc.R script, the function cross.analyze() is used to obtain which trios were significant     from GMAC for each tissue (after Qvalue correction). These trios are then used as candidates for the simulations described   below.

- In each scenario, we simulate the trans gene T2* under a specific model i.e 

		T2* ~ ^b0 + ^b1*V + ^b2*T1 + ^b3*x1 + ^b4*x2 + ... + N(0, ^sigma)

  where each xi is a confouding variable 
  ^ denotes an estimated value (i.e hat)
  and the residuals are simulated from a normal distribution centered at zero with variance == estimated resiual SE
  and the slope coefficients + intercept are the estimated coefficients from fitting each trans gene to the model 
		
		T2 ~ b0 + b1*V + b2*T1 + b3*x1 + b4*x2 + ... + E

- In GMAC_postproc.R the function to execute all simulations is simu12(). It is a wrapper function that wraps simu1(), simu2(), simu3(), and simu4().

- The R-script get_median_p.R is used to obtain the median pvalue from repeated tests under each simulation scenario.

We looked at four possible scenarios 
	(i) When the inference model is overspecificed compared to the true model (true model has fewer confounders)
	    - we simulate the true trans gene using all selected confounders with a pvalue below alpah = 0.001. We then  	              use all selected confounders in the inference model 
	    - In GMAC_postproc.R the function to execute this simulation is simu1()


	(ii) When the inference model is underspecified compared to the true model (true model has more confounders)
	     - we simulate the true trans gene using all selected confounders plus a random number of additional confounders
               (we determine the number of additional WB pcs by selecting from a Uniform[1, 20]). We then analyze using only
               the selected PCs
	     - In GMAC_postproc.R the function to execute this simulation is simu2()


	(iii) When the inference model includes only the confounders from from MRPC and the true model is the 
	      GMAC model (extreme underspecification)
	      - We simulate the true model using the GMAC selected PCs and then infer using the MRPC selected PCs
              - In GMAC_postproc.R the function to execute this simulation is simu3()


	(iv) When the true model and inference model are both the GMAC model (control scenario to understand the stability 		     of permutation)
	     - We simulate the true model with the GMAC selected PCs and then infer the using the same (correct) model. 
             - In GMAC_postproc.R the function to execute this simulation is simu4()


- All simulation scripts are stored in MRPC_support/GMAC Analysis/sim_scripts/

- All simulation master tables are stored in MRPC_support/GMAC Analysis/master_tables/

**get_median_sim_p.R - contains the script for run the 4 sim scenarios and calc the median p-value