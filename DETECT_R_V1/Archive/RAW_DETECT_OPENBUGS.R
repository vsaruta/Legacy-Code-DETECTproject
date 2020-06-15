model
{
  # HERE double check "N" for each iterative
  i_n = 4
  t_n = Nt
  z_n = Nz
  for (i in 1:i_n)
  {
    for (z in 1:z_n)
    {
      # Likelihood for CO2 data; 
      # CO2_Emi[i] ~ dnorm(CO2_Emission[i,z], tau_CO2_Emi[i,z])
      CO2_flux[i,z] ~ dnorm(cctype[i,z], tau_cctype[t])
      # Compute squared deviations for computing posterior predictive loss:
      #squared_deviations[i] <- pow(CO2_Emi[i]-CO2_Emi.rep[i],2)
      for (t in 1:t_n)
      {
        # Numerical approximation of CO2 concentrations based on the eq. 11 from
        # Ryan et al. 2018:
        #Make sure to sum over types (C12 and C13 and add together microbs and roots)
        # HERE Stype_Star should be recalulated to sum up over ...(one of demenations?)
        # HERE check with notes from Kiona based on Matlab Kim code
        Stype_star = 1/100 * (sum(Stype[i,t,z+1]))
        
        # HERE catm may no need in [t] Need to check
        # HERE Make sure units is CM and then translate to M
        cctype[i,z] <- (Stype_star[t]/Dgs_t[t,z])*(z-((z^2)/2))+upbc[t]
        
      }
    }
  }
  for(t in 1:t_n){
    tau_cctype[t] ~ dgamma(0.1, 0.1)
  }
}

