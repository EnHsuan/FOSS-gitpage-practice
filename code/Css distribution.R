library(httk)

chem <- read.csv("19 overlapping chemical list.csv")

set.seed(1234)

DT_pop <- httkpop_generate(method = "d", #direct resampling
                           nsamp = 10000, #number of individuals
                           gfr_resid_var = TRUE,
                           ckd_epi_race_coeff = FALSE)
DT_mc <- list()
css <- data.frame(matrix(NA, nrow=10000, ncol=19))
colnames(css) <- chem$CAS.No.
par.urine <- list()
fub <- data.frame(matrix(NA, nrow=10000, ncol=19))
gfr <- data.frame(matrix(NA, nrow=10000, ncol=19))
bsa <- data.frame(matrix(NA, nrow=10000, ncol=19))
bw <- data.frame(matrix(NA, nrow=10000, ncol=19))
for (i in 1:19){
  cas <- chem[i,58]
  DT_mc[[i]] <- create_mc_samples(chem.cas=cas, species="Human",
                                  model = "3compartmentss", #insert your chosen model here
                                  httkpop.dt = DT_pop,
                                  httkpop = TRUE, #in httk v2.2.2, this will no longer overwrite your user-supplied 
                                                  #httkpop.dt value – and now it's needed to tell create_mc_samples() 
                                                  #that you are using population physiology
                                  samples = nrow(DT_pop) #need to tell create_mc_samples() how many samples you are 
                                                         #giving it, so it can generate the correct number of samples 
                                                         #for other parameters
  )
  css[, i] <- DT_mc[[i]][,
                         calc_analytic_css(chem.cas=cas, output.units='mg/L', species="Human",
                                           parameters = .SD , #.SD = the subset of columns of DT_mc specified in .SDcols argument
                                           model = "3compartmentss",
                                           suppress.messages = TRUE),
                         by = 1:nrow(DT_mc[[i]]), #go row by row
                         .SDcols = names(DT_mc[[i]]) #pass all columns when .SD is used, not just a subset of them
  ]$V1
  
  
  #extract GFR and fraction unbound, also go row by row
  colnames(DT_mc[[i]])[1] <- "weight_adj"
  par.urine[[i]] <- merge(DT_mc[[i]], DT_pop, by="weight_adj")
  fub[, i] <- par.urine[[i]]$Funbound.plasma
  gfr[, i] <- par.urine[[i]]$gfr_est
  bsa[, i] <- par.urine[[i]]$BSA_adj
  bw[, i] <- par.urine[[i]]$weight_adj
}


#Css
colnames(css) <- chem$CAS.y
write.csv(css, file="Css samples.csv", row.names = FALSE)

css.t <- as.data.frame(t(css))
css.t$CAS <- row.names(css.t)
write.csv(css.t, file="Css samples CAS by row.csv", row.names = FALSE)


css.ln <- log(css)
library(reshape2)
css.ln.m <- melt(css.ln)

library(ggplot2)
plot <- ggplot(css.ln.m, aes(sample = value))+
  geom_qq()+
  geom_qq_line(color="red")+
  facet_wrap(~ `variable`, scales = "free")+
  theme_bw()
ggsave(plot, file="mc samples qqplot lognormal transformed.pdf", width=20, height=20, scale=0.5)

css.q <- data.frame(t(apply(css, 2, quantile, probs=c(0.5,0.05,0.95), na.rm=TRUE)))
css.q$CAS <- row.names(css.q)
css.plot <- ggplot(css.q, aes(x=X50., y=CAS))+
  geom_point()+
  geom_errorbar(aes(xmin=X5., xmax=X95.))+
  scale_x_log10()+
  theme_bw()
print(css.plot)
ggsave(css.plot, file="Css quantile plot.pdf", width=8, height=10)

#GFR (original unit: mL/min/1.73 m^2 body surface area)
colnames(gfr) <- chem$CAS.y
write.csv(gfr, file="GFR samples.csv", row.names = FALSE)

gfr.t <- as.data.frame(t(gfr))
gfr.t$CAS <- row.names(gfr.t)
write.csv(gfr.t, file="GFR samples CAS by row.csv", row.names = FALSE)

#Fraction unbound (unitless
colnames(fub) <- chem$CAS.y
write.csv(fub, file="Fub samples.csv", row.names = FALSE)

fub.t <- as.data.frame(t(fub))
fub.t$CAS <- row.names(fub.t)
write.csv(fub.t, file="Fub samples CAS by row.csv", row.names = FALSE)
 
#BSA adjusted (unit: cm^2)
colnames(bsa) <- chem$CAS.y
write.csv(bsa, file="BSA adjusted samples.csv", row.names = FALSE)

bsa.t <- as.data.frame(t(bsa))
bsa.t$CAS <- row.names(bsa.t)
write.csv(bsa.t, file="BSA adjusted samples CAS by row.csv", row.names = FALSE)

#Body weight adjusted (unit: kg)
colnames(bw) <- chem$CAS.y
write.csv(bw, file="Body weight adjusted samples.csv", row.names = FALSE)

bw.t <- as.data.frame(t(bw))
bw.t$CAS <- row.names(bw.t)
write.csv(bw.t, file="Body weight adjusted samples CAS by row.csv", row.names = FALSE)

#Uss_tk (unit: L/kg/day)
uss_tk <- gfr*(1/1000)*1440*(bsa/10000)*fub*(1/bw)
colnames(uss_tk) <- chem$CAS.y
write.csv(uss_tk, file="GFRxFub samples.csv", row.names = FALSE)

uss_tk.t <- as.data.frame(t(uss_tk))
uss_tk.t$CAS <- row.names(uss_tk.t)
write.csv(uss_tk.t, file="GFRxFub samples CAS by row.csv", row.names = FALSE)
