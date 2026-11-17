
library(dplyr)
library(sf)
library(terra)
library(tidyverse)
#library(factoextra)
#library(lmtest)
#library(tidyr)
library(mgcv)
library(gridExtra)
#library(marginaleffects)
library(DHARMa)
library(ggplot2)
library(GGally)
library(rnaturalearth)

setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder/")

atmoAnnotation_tom <- read.csv("outputData/ForGiulia_backgandpresence_forAnnotation_Jul2025_annotated.csv")

# Load model dataset (sf object)
model_data <- readRDS("outputData/PCA_wind_pred_2708cio.rds")

head(model_data)
names(model_data)

# store laea coordinates
crsLaea <- crs(model_data)
model_data$x_laea <- st_coordinates(model_data)[,1]
model_data$y_laea <- st_coordinates(model_data)[,2]

# Reproject to swiss coordinate system LV95 (EPSG:2056)
model_data <- st_transform(model_data, 2056)
crsLv95 <- crs(model_data)
model_data$x_lv95 <- st_coordinates(model_data)[,1]
model_data$y_lv95 <- st_coordinates(model_data)[,2]

# extent to download swiss dem
st_bbox(model_data)

# map data
switz <- ne_countries(scale = "medium", country = "Switzerland", returnclass = "sf")
ggplot() +
  geom_sf(data=switz) +
  geom_sf(data=model_data, aes(col=territory_ID, alpha=as.factor(response))) +
  theme_void()

# number of nests per territory
nestsPerTerr <- model_data %>%
  group_by(territory_ID) %>%
  summarise(Npres = sum(response == 1), Nabs  = sum(response == 0), 
            ratio_abs_pres = Nabs / Npres) %>% 
  st_drop_geometry()
print(nestsPerTerr, n=Inf)

# correlation values between interesting variables (what is the updraft coefficient mean??)
# means and variances separately
vars_mean <- c("slope", "dtm","ASHFL_S_mean", "w_oro_mean", "w_oro_uplift_30m_mean", "updraft_coeff_mean", 
          "U_mean", "V_mean", "W_mean", "aspect", "TKE_mean", "distance_edge_territory")
varsComb <- expand.grid(var1 = vars_mean, var2 = vars_mean, stringsAsFactors = FALSE)

# remove self-pairs and keep only one of each unordered pair
varsComb <- varsComb %>%
  filter(var1 != var2)
df <- st_drop_geometry(model_data)[, vars_mean]
cor_result <- apply(varsComb, 1, function(row){
  cor(df[[row["var1"]]], df[[row["var2"]]],
      use = "pairwise.complete.obs",
      method = "pearson")
})
cor_table <- cbind(varsComb, cor = cor_result)

# Check high correlations. We exclude updraft_coeff_mean and will have to decide between W_mean and w_oro_mean
cor_table[cor_table$cor>0.4,]
# Check correlation between the three w_oro. Why the two oro uplift aren't correlated?
filter(cor_table, var1=="w_oro_mean" & var2=="updraft_coeff_mean")
filter(cor_table, var1=="w_oro_uplift_30m_mean" & var2=="updraft_coeff_mean")
filter(cor_table, var1=="w_oro_mean" & var2=="w_oro_uplift_30m_mean")

# Same for other summaries, variances and maximum values
#vars_var <- c("ASHFL_S_variance", "w_oro_variance", "w_oro_uplift_30m_mean", "U_variance", "V_variance", "W_variance", "TKE_variance")
vars_max <- c("ASHFL_S_max","w_oro_max","w_oro_uplift_30m_mean","U_max","V_max","W_max","TKE_max")
vars <- c(vars_mean, vars_max)

varsComb2 <- expand.grid(var1 = vars, var2 = vars, stringsAsFactors = FALSE)
varsComb2 <- varsComb2 %>%
  filter(var1 != var2)
df <- st_drop_geometry(model_data)[, vars]
cor_result <- apply(varsComb2, 1, function(row){
  cor(df[[row["var1"]]], df[[row["var2"]]],
      use = "pairwise.complete.obs",
      method = "pearson")
})
cor_table2 <- cbind(varsComb2, cor = cor_result)
cor_table2 <- cor_table2[!duplicated(cor_table2),]
# Check high correlations. We will have to decide between W_max and w_oro_max
cor_table2[cor_table2$cor<0.4,]
uncor_vars <- unique(cor_table2$var1[cor_table2$cor<0.4])
uncor_vars2 <- unique(cor_table2$var2[cor_table2$cor<0.4])
unique(c(uncor_vars2, uncor_vars))

# Let's check both w_oro_mean (calculated with 1 km dem) and "w_oro_uplift_30m_mean" (with 30 m dem) just as sanity check
# 
ggpairs(
  data = st_drop_geometry(model_data[ , vars_mean]),
  upper = list(continuous = wrap("points", alpha = 0.6, size = 0.5)),
  lower = list(continuous = wrap("cor", size = 5)),
  diag = list(continuous = wrap("densityDiag"))
)
ggpairs(
  data = st_drop_geometry(model_data[ , vars_var]),
  upper = list(continuous = wrap("points", alpha = 0.6, size = 0.5)),
  lower = list(continuous = wrap("cor", size = 5)),
  diag = list(continuous = wrap("densityDiag"))
)

# Check average slope and dem of nest locations
hist(model_data$dtm[model_data$response==1])
hist(model_data$slope[model_data$response==1])
plot(slope~dtm, data=model_data[model_data$response==1,])


# Boxplot of the above variables between presences and absences
df_long <- model_data %>%
  select(response, all_of(vars)) %>%
  pivot_longer(
    cols = all_of(vars),
    names_to = "variable",
    values_to = "value")

ggplot(df_long, aes(x = factor(response), y = value, fill = factor(response))) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    x = "Response (0/1)",
    y = NULL,
    fill = "Response",
    title = "Boxplots of predictor variables by response group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "bottom"
  )

# Scale variables
model_data$dtm_std <- scale(model_data$dtm)
model_data$ASHFL_S_mean_std <- scale(model_data$ASHFL_S_mean)
model_data$TKE_mean_std <- scale(model_data$TKE_mean)
model_data$UVspeed_10M_mean_std <- scale(model_data$UVspeed_10M_mean)
model_data$UVdir_10m_mean_std <- scale(model_data$UVdir_10m_mean)
model_data$updraft_coeff_mean_std <- scale(model_data$updraft_coeff_mean)
model_data$w_oro_mean_std <- scale(model_data$w_oro_mean)
model_data$w_oro_uplift_30m_mean_std <- scale(model_data$w_oro_uplift_30m_mean)
model_data$slope_std <- scale(model_data$slope)
model_data$aspect_std <- scale(model_data$aspect)
model_data$distance_edge_territory_std <- scale(model_data$distance_edge_territory)
model_data$TRI_std <- scale(model_data$TRI)
model_data$UVspeed_10M_mean_std <- scale(model_data$UVspeed_10M_mean)
model_data$U_mean_std <- scale(model_data$U_mean)
model_data$V_mean_std <- scale(model_data$V_mean)
model_data$W_mean_std <- scale(model_data$W_mean)

# Old method, using gam with cloglog and weights, only one model with all variables
# Adjust weights: 
# Much more weights to absences than to presences to model presence only data (Renner et al MEE 2015)
# model_data$cellweights <- ifelse(model_data$response == 1, 1, 1E4)
# best_mod_w <- mgcv::gam(response ~ 
#                           s(UVdir_10m_mean_std, bs="cc", k=16)
#                         + s(UVspeed_10M_mean_std, bs="tp", k=8)
#                         + s(slope_std, bs="tp", k=9)
#                         + s(aspect_std, bs="cc", k=9)
#                         + s(ASHFL_S_mean_std, bs="tp", k=19)
#                         + s(w_oro_uplift_30m_mean_std, bs="tp", k=7)
#                         + s(TKE_mean_std, bs="tp", k=7)
#                         + s(dtm_std, bs="tp", k=6)
#                         + ti(UVdir_10m_mean_std, UVspeed_10M_mean_std, bs=c("cc", "tp"), k= c(16, 8))
#                         + s(territory_ID, bs="re")
#                         + s(x, y, bs="tp", k=120),
#                         weights = cellweights,
#                         data = model_data,
#                         family = binomial(link = "cloglog"), method = "fREML", discrete = T)
# Bjorn says that this is not a real conditional logistic regression, we better use the method below suggested by Natasha Klappstein


# New method:
# Multiple-case matched sets (multiple presence and absence per stratum)
# We use conditional logistic regression following Klappstein et al ...
# ~ rule of thumb 1 predictor per 10 nests
table(model_data$response)
str(model_data)

## Make a first model with only vertical distance, distance from territory border, slope, dem and land cover
## predict suitable cells within each territory
## subset absences to only suitable areas
## run second model with only aspect and atmo variables

time <- 1
best_mod_clogit <- mgcv::gam(cbind(time, as.factor(territory_ID)) ~ 
                          slope_std
                        + s(aspect_std, bs="cc")
                        + V_mean_std
                        + ASHFL_S_mean_std
                        + w_oro_uplift_30m_mean_std
                        + TKE_mean_std,
                        # + s(dtm_std, bs="tp", k=6)
                        # + s(U_mean_std, V_mean_std),
                        weights = response,
                        data = model_data,
                        family = cox.ph)
best_mod_clogit <- mgcv::gam(cbind(time, as.factor(territory_ID)) ~ 
                               slope_std +
                               distance_edge_territory_std,
                             weights = response,
                             data = model_data,
                             family = cox.ph)
best_mod_clogit <- mgcv::gam(cbind(time, as.factor(territory_ID)) ~ 
                            slope_std +
                             dtm_std,
                             weights = response,
                             data = model_data,
                             family = cox.ph)
summary(best_mod_clogit)
plot(best_mod_clogit)
AIC(best_mod_clogit)
draw(best_mod_w)
resid <- simulateResiduals(fittedModel = best_mod_w, plot=F)
plot(resid)








library(survival)
library(mgcv)

sc <- mgcv::smoothCon(te(slope, dtm), data = model_data)[[1]]
X_te <- as.data.frame(sc$X)
colnames(X_te) <- sc$label

dat2 <- cbind(model_data, X_te)
best_mod_clogit <- coxph(
  Surv(rep(1, nrow(dat2)), response) ~ 
    . + strata(territory_ID),
  data = dat2[, c("response", "territory_ID", colnames(X_te))]
)



devtools::install_github("Schaubert/ClogitTree")
library(CLogitTree)
data(illu.small)
head(illu.small)

# Compute CLogitForest based on 4 trees
set.seed(1860)
illu.rf <- CLogitForest(illu.small, response = "y", 
                        exposure = "x", s = "strata", 
                        ntree = 4, depth_max = 2, 
                        tune.mtry = FALSE)

# show estimated exposure effect
illu.rf
p <- predict(illu.rf, illu.small)

summary(illu.rf)

# compute and plot variable importance
(vi <- varimp(illu.rf))
plot(vi)







