###############################################################################
# Multi-species occupancy model fitted to full dataset (2017–2023)
# Four-dimensional structure: Site × Year × Session × Species
# Sessions are defined within each camera-year
# Analysis restricted to the first two sessions per camera-year
# Covariates: DistWater, DistVillag, TreeCover, DogFreq, HumanFreq,
#             SmallMammals, Patch Size
# Categorical predictor: Habitat
###############################################################################

setwd("C:/...")

library(jagsUI)
library(dplyr)
library(tidyr)

###############################################################################
# 1) LOAD DATA AND DEFINE SAMPLING SESSIONS
###############################################################################

DF <- read.csv("BData.csv", header=TRUE)

# Sort observations to ensure consistent indexing
DF <- DF[order(DF$Camera, DF$Year, DF$Survey), ]

# Define temporal sampling sessions nested within each camera-year combination
DF$Session <- ave(
  DF$Camera,
  DF$Camera, DF$Year,
  FUN = seq_along
)

# Restrict analysis to the first two sampling sessions per camera-year
DF <- DF %>% filter(Session %in% c(1,2))

str(DF)
table(DF$Session)

###############################################################################
# 2) DEFINE LEVELS AND CORE DIMENSIONS
###############################################################################

site_levels  <- sort(unique(DF$Camera))
year_levels  <- sort(unique(DF$Year))
nsites       <- length(site_levels)
nyears       <- length(year_levels)
nsess        <- 2

species_cols <- which(names(DF) %in% c(
  "Fishing.Cat","Leopard.Cat","Smooth.coated.otter","Hairy.nosed.otter",
  "Large.spotted.civet","Common.Palm.civet",
  "Long.tailed.macaque","Sunda.Pangolin"
))
(species_names <- names(DF)[species_cols])
nspecies <- length(species_cols)

###############################################################################
# 3) CONSTRUCT FOUR-DIMENSIONAL ARRAYS
#    Detection data (y_4D) and sampling effort (nOcc_4D)
###############################################################################

y_4D <- array(0L,
              dim=c(nsites, nyears, nsess, nspecies),
              dimnames=list(Site=site_levels,
                            Year=year_levels,
                            Sess=1:nsess,
                            Species=species_names))

nOcc_4D <- array(0L,
                 dim=c(nsites, nyears, nsess, nspecies),
                 dimnames=list(Site=site_levels,
                               Year=year_levels,
                               Sess=1:nsess,
                               Species=species_names))

# Populate detection and effort arrays using camera, year, and session indices
for(i in 1:nrow(DF)){
  si <- match(DF$Camera[i], site_levels)
  ti <- match(DF$Year[i],   year_levels)
  ui <- DF$Session[i]

  y_4D[si,ti,ui,]    <- as.integer(DF[i, species_cols])
  nOcc_4D[si,ti,ui,] <- DF$nOcc[i]
}

###############################################################################
# 4) STACKING ACROSS DIMENSIONS
#    Only site–year–session–species combinations with positive sampling effort
###############################################################################

idx <- which(nOcc_4D > 0, arr.ind=TRUE)

site_idx <- idx[,1]
year_idx <- idx[,2]
sess_idx <- idx[,3]
sp_idx   <- idx[,4]

Nobs <- nrow(idx)

y_stack    <- y_4D[idx]
nOcc_stack <- nOcc_4D[idx]

###############################################################################
# 5) SITE × YEAR CONTINUOUS COVARIATES
#    Log-transformed (where appropriate) and standardized
###############################################################################

scale_num <- function(x) as.numeric(scale(x)[,1])

DF$DistWater_sc  <- scale_num(log1p(DF$DistWater))
DF$DistVill_sc   <- scale_num(log1p(DF$DistVillag))
DF$TreeCov_sc    <- scale_num(DF$Tree.cover)
DF$NDVI_sc       <- scale_num(DF$NDVI)
DF$Dog_sc        <- scale_num(DF$Dog.Freq)
DF$Human_sc      <- scale_num(DF$Human.Freq)
DF$SmallM_sc     <- scale_num(DF$Freq.small.mammals)
DF$Patch_sc      <- scale_num(DF$Patch)

# Helper function to build site × year matrices
build_mat <- function(df, col){
  M <- matrix(NA_real_, nsites, nyears,
              dimnames=list(site_levels, year_levels))
  for(i in 1:nrow(df)){
    si <- match(df$Camera[i], site_levels)
    ti <- match(df$Year[i],   year_levels)
    M[si,ti] <- df[[col]][i]
  }
  M[is.na(M)] <- 0
  M
}

distwater_mat <- build_mat(DF, "DistWater_sc")
distvill_mat  <- build_mat(DF, "DistVill_sc")
treecover_mat <- build_mat(DF, "TreeCov_sc")
ndvi_mat      <- build_mat(DF, "NDVI_sc")
dog_mat       <- build_mat(DF, "Dog_sc")
human_mat     <- build_mat(DF, "Human_sc")
smallm_mat    <- build_mat(DF, "SmallM_sc")
patch_mat     <- build_mat(DF, "Patch_sc")

# Extract stacked covariate values matching stacked observations
stack_sy <- function(mat) mat[cbind(site_idx, year_idx)]

distwater_stack <- stack_sy(distwater_mat)
distvill_stack  <- stack_sy(distvill_mat)
treecover_stack <- stack_sy(treecover_mat)
ndvi_stack      <- stack_sy(ndvi_mat)
dog_stack       <- stack_sy(dog_mat)
human_stack     <- stack_sy(human_mat)
smallm_stack    <- stack_sy(smallm_mat)
patch_stack     <- stack_sy(patch_mat)

###############################################################################
# 6) HABITAT (CATEGORICAL SITE × YEAR PREDICTOR)
###############################################################################

hab_levels <- sort(unique(DF$Habitat))

hab_mat <- matrix(NA_integer_, nsites, nyears,
                  dimnames=list(site_levels, year_levels))

for(i in 1:nrow(DF)){
  si <- match(DF$Camera[i], site_levels)
  ti <- match(DF$Year[i],   year_levels)
  hab_mat[si,ti] <- match(DF$Habitat[i], hab_levels)
}

hab_stack <- hab_mat[cbind(site_idx, year_idx)]
nHab <- length(hab_levels)

###############################################################################
# JAGS MODEL DEFINITION
###############################################################################

cat(file = "stacked_4D.txt", "
model{

  #######################################################
  # 1. LIKELIHOOD: STACKED OBSERVATIONS
  #    (implicit 4D indexing)
  #######################################################

  for(s in 1:Nobs){

    # Occupancy model (logit scale)
    lpsi_s[s] <-
        mu.lpsi +
        eps_psi[ sp_idx[s] ] +
        a_year[ year_idx[s] ] +
        gamma_sp[ sp_idx[s] ] * year_c[ year_idx[s] ] +

        # Site × year covariates, replicated across sessions
        beta_distwater[ sp_idx[s] ] * distwater[s] +
        beta_distvill[ sp_idx[s] ] * distvill[s] +
        b_treecover * treecover[s] +

        # Habitat-specific fixed effects
        alpha_hab[ hab_idx[s] ]

    psi_s[s] <- ilogit(lpsi_s[s])
    z[s] ~ dbern(psi_s[s])

    # Detection model (logit scale)
    lp_s[s] <-
        mu.lp +
        eps_p[ sp_idx[s] ] +
        b_year[ year_idx[s] ]

    p_s[s] <- ilogit(lp_s[s])

    # Binomial observation model conditional on occupancy
    y[s] ~ dbin( p_s[s] * z[s] , nrep[s] )

  }

  #######################################################
  # 2. RANDOM EFFECTS
  #######################################################

  for(k in 1:nSpecies){
    eps_psi[k] ~ dnorm(0, tau_sp_psi)
    eps_p[k]   ~ dnorm(0, tau_sp_p)
    gamma_sp[k] ~ dnorm(0, tau_gamma)
  }

  for(t in 1:nYears){
    a_year[t] ~ dnorm(0, tau_year_psi)
    b_year[t] ~ dnorm(0, tau_year_p)
  }

  #######################################################
  # 3. PRIORS – OCCUPANCY
  #######################################################

  psi_mean ~ dbeta(1,1)
  mu.lpsi <- logit(psi_mean)

  b_treecover ~ dlogis(0,1)
  b_distwater ~ dlogis(0,1)
  b_distvill  ~ dlogis(0,1)

  for(k in 1:nSpecies){
    beta_distwater[k] ~ dnorm(b_distwater, tau_beta_distwater)
    beta_distvill[k]  ~ dnorm(b_distvill,  tau_beta_distvill)
  }

  sd_beta_distwater ~ dunif(0,5)
  tau_beta_distwater <- 1 / (sd_beta_distwater*sd_beta_distwater)

  sd_beta_distvill ~ dunif(0,5)
  tau_beta_distvill <- 1 / (sd_beta_distvill*sd_beta_distvill)

  alpha_hab[1] <- 0
  for(h in 2:nHab){
    alpha_hab[h] ~ dnorm(0, tau_hab)
  }

  #######################################################
  # 4. PRIORS – DETECTION
  #######################################################

  p_mean ~ dbeta(1,1)
  mu.lp <- logit(p_mean)

  #######################################################
  # 5. HYPERPRIORS
  #######################################################

  sd_sp_psi ~ dunif(0,5)
  tau_sp_psi <- 1 / (sd_sp_psi*sd_sp_psi)

  sd_sp_p ~ dunif(0,5)
  tau_sp_p <- 1 / (sd_sp_p*sd_sp_p)

  sd_year_psi ~ dunif(0,5)
  tau_year_psi <- 1 / (sd_year_psi*sd_year_psi)

  sd_year_p ~ dunif(0,5)
  tau_year_p <- 1 / (sd_year_p*sd_year_p)

  sd_gamma ~ dunif(0,5)
  tau_gamma <- 1 / (sd_gamma*sd_gamma)

  sd_hab ~ dunif(0,5)
  tau_hab <- 1 / (sd_hab*sd_hab)

}
")

###############################################################################
# 7) YEAR CENTERING
###############################################################################

years  <- 1:nyears
year_c <- years - mean(years)

###############################################################################
# 8) JAGS DATA LIST
###############################################################################

str(jags_data <- list(
  Nobs      = Nobs,
  y         = y_stack,
  nrep      = nOcc_stack,
  site_idx  = site_idx,
  year_idx  = year_idx,
  sess_idx  = sess_idx,
  sp_idx    = sp_idx,
  nSites    = nsites,
  nYears    = nyears,
  nSess     = nsess,
  nSpecies  = nspecies,
  distwater = distwater_stack,
  distvill  = distvill_stack,
  treecover = treecover_stack,
  hab_idx   = hab_stack,
  nHab      = nHab,
  year_c    = year_c
))

###############################################################################
# 9) INITIAL VALUES
###############################################################################

inits <- function() list(
  z = ifelse(jags_data$y>0,1,rbinom(Nobs,1,0.3)),
  psi_mean=0.5,
  p_mean=0.5,
  sd_sp_psi  = runif(1,0.1,1.5),
  sd_sp_p    = runif(1,0.1,1.5),
  sd_year_psi= runif(1,0.1,1.5),
  sd_year_p  = runif(1,0.1,1.5),
  sd_gamma   = runif(1,0.1,1.5),
  sd_hab     = runif(1,0.1,2)
)

###############################################################################
# 10) PARAMETERS MONITORED
###############################################################################

params <- c(
  "psi_mean","p_mean","mu.lpsi","mu.lp",
  "eps_psi","eps_p","gamma_sp",
  "a_year","b_year",
  "sd_sp_psi","sd_sp_p",
  "sd_year_psi","sd_year_p",
  "sd_gamma","sd_hab",
  "b_distwater","b_distvill","b_treecover",
  "beta_distwater","sd_beta_distwater",
  "beta_distvill","sd_beta_distvill",
  "alpha_hab","psi_s","p_s"
)

###############################################################################
# 11) MODEL FITTING
###############################################################################

fit_4D <- jags(
  data=jags_data,
  inits=inits,
  parameters.to.save=params,
  model.file="stacked_4D.txt",
  n.chains=3,
  n.adapt=500000,
  n.iter=1000000,
  n.burnin=500000,
  n.thin=10,
  parallel=TRUE
)

print(fit_4D, digits=3)
traceplot(fit_4D)

###############################################################################
# 12) EXTRACTION OF POSTERIOR OCCUPANCY ESTIMATES
###############################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

m <- as.matrix(fit_4D$samples)

psi_cols <- grep("^psi_s\\[[0-9]+\\]$", colnames(m), value=TRUE)
psi_m <- m[, psi_cols, drop=FALSE]

iters_keep <- sample(1:nrow(psi_m), size = 1000)
psi_m <- psi_m[iters_keep, , drop=FALSE]

idx_s <- as.integer(str_match(psi_cols, "psi_s\\[([0-9]+)\\]")[,2])

psi_all <- psi_m %>%
  as.data.frame() %>%
  mutate(iter=row_number()) %>%
  pivot_longer(cols=all_of(psi_cols),
               names_to="psi_name", values_to="psi")

###############################################################################
# 13) VISUALIZATION OF OCCUPANCY BY SPECIES AND HABITAT
###############################################################################


psi_clean <- psi_all %>%
    mutate(
        # limpiar nombres: de "Fishing.cat" → "Fishing Cat"
        species_label_clean = gsub("\\.", " ", species_label),
        species_label_clean = str_to_title(species_label_clean),

        habitat_label_clean = habitat_label   # ya viene limpio
    ) %>%
    filter(!is.na(species_label_clean),
           !is.na(habitat_label))


set.seed(123)

psi_fast <- psi_clean %>%
  group_by(habitat_label_clean, species_label_clean) %>%
  sample_frac(0.20) %>%
  ungroup()


species_colors_new <- c(
  "Fishing Cat"         = "#E3C18F",
  "Leopard Cat"         = "#C19A1B",
  "Smooth Coated Otter" = "#5A3E2B",
  "Hairy Nosed Otter"   = "#5F9EA0",
  "Large Spotted Civet" = "#B8860B",
  "Common Palm Civet"   = "#6B8E23",
  "Long Tailed Macaque" = "#A0522D",
  "Sunda Pangolin"      = "#708090"
)

psi_fast$species_label_clean <- factor(
    psi_fast$species_label_clean,
    levels = names(species_colors_new)
)

psi_fast$habitat_label_clean <- factor(
    psi_fast$habitat_label_clean,
    levels = sort(unique(psi_fast$habitat_label_clean))
)


dev.new(width=9, height=9)

g <- ggplot(
        psi_fast,
        aes(x = species_label_clean, y = psi, fill = species_label_clean)
    ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.88, color="grey35") +
    stat_summary(fun = median, geom="point",
                 shape=23, fill="white", size=1.6) +
    facet_wrap(~ habitat_label_clean, ncol=2, scales="free_y") +
    labs(title = "", x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
        plot.margin = margin(0.05, 0.4, 0.4, 0.4, "cm"),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = c(0.86, 0.04),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 14),
        legend.box = "horizontal",
        legend.key.width = unit(18, "pt")) +
    scale_fill_manual(
      values = species_colors_new,
      limits = names(species_colors_new)  # <- fija orden de la leyenda
    )+
    guides(fill = guide_legend(title = NULL))

print(g)
