holder <- strategies_of_interest[! strategies_of_interest %in% c("0.9999 0 HLA plus IgATTG" , 
                                                                 "0.882 0.417 HLA plus IgATTG",
                                                                 "0.807 0.61 HLA plus IgATTG" ,
                                                                 "0.667 0.872 HLA plus IgATTG",
                                                                 "0.533 0.952 HLA plus IgATTG" ,
                                                                 "0.331 0.987 HLA plus IgATTG")]

m_of_interest <- bcea(e = t(output$total_qalys[holder,]), 
                      c = t(output$total_costs[holder,]), ref = 1, 
                      interventions = holder)
summary(m_of_interest) 



mce <- multi.ce(m_of_interest)

ceac.plot(mce, graph = c("base"), pos = "topright",
          cex.main = 0.2,
          line = list(colors = c(1:length(holder)))) #CEAC 

