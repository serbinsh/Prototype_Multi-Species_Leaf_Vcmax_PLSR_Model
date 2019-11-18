####################################################################################################
#
#  
#   Predict Leaf Vcmax using a prototype multi-species PLSR model (Serbin et al. unpublished) applied 
#   to fresh leaf spectra collected on aspen and cottonweed seedlings grown in the UW-Madison Biotron 
#   facility 
#
#   Spectra and trait data source:
#   https://ecosis.org/package/2008-university-of-wisconsin-biotron-fresh-leaf-spectra-and-gas-exchange-leaf-traits
#
#    Notes:
#    * Provided as a basic example of how to apply the model to new spectra observations
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
#
#    --- Last updated:  11.14.2019 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

list.of.packages <- c("readr","scales","plotrix","httr","devtools")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load libraries needed for script
library(readr)    # readr - read_csv function to pull data from EcoSIS
library(plotrix)  # plotCI - to generate obsvered vs predicted plot with CIs
library(scales)   # alpha() - for applying a transparency to data points
library(devtools)
library(httr)

# define function to grab PLSR model from GitHub
#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <-function(url, sep = ",", header = TRUE) {
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory (scratch space)
output_dir <- file.path("~",'scratch/')
if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE, showWarnings = FALSE)
setwd(output_dir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### PLSR Coefficients - Grab from GitHub
print("**** Downloading PLSR coefficients ****")
git_repo <- "https://raw.githubusercontent.com/serbinsh/Prototype_Multi-Species_Leaf_Vcmax_PLSR_Model/master/"
githubURL <- paste0(git_repo,"PLSR_model_coefficients/Vcmax_PLSR_Coefficients_9comp.csv")
LeafVcmax.plsr.coeffs <- source_GitHubData(githubURL)
rm(githubURL)
githubURL <- paste0(git_repo,"PLSR_model_coefficients/Leaf_Vcmax_Jackkife_PLSR_Coefficients.csv")
LeafVcmax.plsr.jk.coeffs <- source_GitHubData(githubURL)
rm(githubURL)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Example datasets
# 
# URL:  https://ecosis.org/package/2008-university-of-wisconsin-biotron-fresh-leaf-spectra-and-gas-exchange-leaf-traits
#
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Grab data
print("**** Downloading Ecosis data ****")
ecosis_id <- "03e46f54-7d68-4a8c-896a-fcea87e9cf10"  # Biotron data
ecosis_file <- sprintf(
  "https://ecosis.org/api/package/%s/export?metadata=true",
  ecosis_id
)
message("Downloading data...")
dat_raw <- read_csv(ecosis_file)
message("Download complete!")
head(dat_raw)
names(dat_raw)[1:40]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Create validation dataset
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)

spectra <- dat_raw[,names(dat_raw)[match(seq(Start.wave,End.wave,1),names(dat_raw))]]
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot data
waves <- seq(500,2400,1)
cexaxis <- 1.5
cexlab <- 1.8
ylim <- 65
ylim2 <- 65

mean_spec <- colMeans(spectra[,which(names(spectra) %in% seq(Start.wave,End.wave,1))])
spectra_quantiles <- apply(spectra[,which(names(spectra) %in% seq(Start.wave,End.wave,1))],
                           2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))

print("**** Plotting Ecosis data. Writing to scratch space ****")
png(file=file.path(output_dir,'UW-Biotron_leaf_spectra_summary_plot.png'),height=3000,
    width=3900, res=340)
par(mfrow=c(1,1), mar=c(4.5,5.7,0.3,0.4), oma=c(0.3,0.9,0.3,0.1)) # B, L, T, R
plot(waves,mean_spec*100,ylim=c(0,ylim),cex=0.00001, col="white",xlab="Wavelength (nm)",
     ylab="Reflectance (%)",cex.axis=cexaxis, cex.lab=cexlab)
polygon(c(waves ,rev(waves)),c(spectra_quantiles[6,]*100, rev(spectra_quantiles[2,]*100)),
        col="#99CC99",border=NA)
lines(waves,mean_spec*100,lwd=3, lty=1, col="black")
lines(waves,spectra_quantiles[1,]*100,lwd=1.85, lty=3, col="grey40")
lines(waves,spectra_quantiles[7,]*100,lwd=1.85, lty=3, col="grey40")
legend("topright",legend=c("Mean reflectance","Min/Max", "95% CI"),lty=c(1,3,1),
       lwd=c(3,3,15),col=c("black","grey40","#99CC99"),bty="n", cex=1.7)
box(lwd=2.2)
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Applying PLSR model to estimate Vcmax from spectral observations ****")
# setup model
dims <- dim(LeafVcmax.plsr.coeffs)
LeafVcmax.plsr.intercept <- LeafVcmax.plsr.coeffs[1,]
LeafVcmax.plsr.coeffs <- data.frame(LeafVcmax.plsr.coeffs[2:dims[1],])
names(LeafVcmax.plsr.coeffs) <- c("wavelength","coefs")
LeafVcmax.plsr.coeffs.vec <- as.vector(LeafVcmax.plsr.coeffs[,2])

# estimate Vcmax
Start.wave <- 500
End.wave <- 2400
sub_spec <- as.matrix(droplevels(spectra[,which(names(spectra) %in% seq(Start.wave,End.wave,1))]))
temp <- as.matrix(sub_spec) %*% LeafVcmax.plsr.coeffs.vec
LeafVcmax <- data.frame(rowSums(temp))+LeafVcmax.plsr.intercept[,2]
LeafVcmax <- LeafVcmax[,1]
LeafVcmax[LeafVcmax<0] <- NA
names(LeafVcmax) <- "FS_PLSR_LeafVcmax_umol_m2_s"
hist(LeafVcmax)

# organize output
LeafVcmax.PLSR.dataset <- data.frame(sample_info, FS_PLSR_LeafVcmax_umol_m2_s=LeafVcmax)
head(LeafVcmax.PLSR.dataset)

# Derive PLSR Vcmax estimate uncertainties - ugly, we can do better than this
print("**** Deriving uncertainty estimates ****")
dims <- dim(LeafVcmax.plsr.jk.coeffs)
intercepts <- LeafVcmax.plsr.jk.coeffs[,2]
jk.leaf.vcmax.est <- array(data=NA,dim=c(dim(sub_spec)[1],dims[1]))
for (i in 1:length(intercepts)){
  coefs <- unlist(as.vector(LeafVcmax.plsr.jk.coeffs[i,3:dims[2]]))
  temp <- sub_spec %*% coefs
  values <- data.frame(rowSums(temp))+intercepts[i]
  jk.leaf.vcmax.est[,i] <- values[,1]
  rm(temp)
}

jk.leaf.vcmax.est.quant <- apply(jk.leaf.vcmax.est,1,quantile,probs=c(0.025,0.975))
jk.leaf.vcmax.est.quant2 <- data.frame(t(jk.leaf.vcmax.est.quant))
names(jk.leaf.vcmax.est.quant2) <- c("FS_PLSR_LeafVcmax_L5","FS_PLSR_LeafVcmax_U95")
jk.leaf.vcmax.est.sd <- apply(jk.leaf.vcmax.est,1,sd)
names(jk.leaf.vcmax.est.sd) <- "FS_PLSR_LeafVcmax_umol_m2_s_Sdev"

## Combine into final dataset
stats <- data.frame(jk.leaf.vcmax.est.sd,jk.leaf.vcmax.est.quant2)
names(stats) <- c("FS_PLSR_LeafVcmax_umol_m2_s_Sdev","FS_PLSR_LeafVcmax_L5","FS_PLSR_LeafVcmax_U95")
LeafVcmax.PLSR.dataset.out <- data.frame(LeafVcmax.PLSR.dataset,stats,
                                     residual=(LeafVcmax.PLSR.dataset$FS_PLSR_LeafVcmax_umol_m2_s - 
                                                 LeafVcmax.PLSR.dataset$Vcmax_PLSR_Temperature))
head(LeafVcmax.PLSR.dataset.out)

# output results
write.csv(x = LeafVcmax.PLSR.dataset.out, file = file.path(output_dir,"UW-Biotron_PLSR_estimated_Vcmax_data.csv"),
          row.names = F)
# calculate error stats
rmse <- sqrt(mean(LeafVcmax.PLSR.dataset.out$residual^2, na.rm=T))
# calculate fit stats
reg <- lm(LeafVcmax.PLSR.dataset.out$FS_PLSR_LeafVcmax_umol_m2_s~LeafVcmax.PLSR.dataset.out$Vcmax_PLSR_Temperature)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot up results
ptcex <- 1.8
cexaxis <- 1.3
cexlab <- 1.8
print("**** Plotting PLSR estimated Vcmax validation plot. Writing to scratch space ****")
png(file=file.path(output_dir,'UW-Biotron_PLSR_estimated_Vcmax_validation_plot.png'),height=3000,
    width=3900, res=340)
par(mfrow=c(1,1), mar=c(4.5,5.4,1,1), oma=c(0.3,0.9,0.3,0.1)) # B, L, T, R
plotCI(LeafVcmax.PLSR.dataset.out$FS_PLSR_LeafVcmax_umol_m2_s,LeafVcmax.PLSR.dataset.out$Vcmax_PLSR_Temperature,
       li=LeafVcmax.PLSR.dataset.out$FS_PLSR_LeafVcmax_L5,gap=0.009,sfrac=0.004,lwd=1.6,
       ui=LeafVcmax.PLSR.dataset.out$FS_PLSR_LeafVcmax_U95,err="x",pch=21,col="black",
       pt.bg=alpha("grey70",0.7),scol="grey30",xlim=c(0,200),cex=ptcex,
       ylim=c(0,200),xlab="Predicted Vcmax (umols/m2/s)",
       ylab="Observed Vcmax (umols/m2/s)",main="",
       cex.axis=cexaxis,cex.lab=cexlab)
abline(0,1,lty=2,lw=2)
legend("topleft",legend = c(paste0("RMSE = ",round(rmse,2)),
                            paste0("R2 = ",round(summary(reg)$r.squared,2))), bty="n", cex=1.5)
box(lwd=2.2)
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
### EOF