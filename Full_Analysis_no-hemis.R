#Final Analysis Pipeline for all analyses of granularity project
library(tidyverse)
library(emmeans)
library(afex)
library(lme4)
library(PupillometryR)
library(merTools)


theme_set(theme_minimal(base_size = 20))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(),
             legend.position = "none",
             axis.title.x = element_blank(),
             axis.title.y = element_text(size = 24))
          
acq.names = list(
  '1400'='Acquisition 1','CAP'='Acquisition 2')
hemi_labeller <- function(variable,value){
  return(acq.names[value])
}

hemi.names = list(
  'r'='right','l'='left')
hemisphere_labeller <- function(variable,value){
  return(hemi.names[value])
}


dodgehalf <- position_dodge(width=0.5)
dodgewhole <- position_dodge(width = 1)

ant.col = viridis::viridis_pal(direction = -1)(4)[1]
ant.col = "#f4c832"
mid.col = viridis::viridis_pal(direction = -1)(4)[2]
mid.col = "#f44f56"
mid_post.col = viridis::viridis_pal(direction = -1)(4)[3]
mid_post.col = "#a270d0"
post.col = viridis::viridis_pal(direction = -1)(4)[4]
post.col = "#3abaee"

ant.mid.col = "#f29b5c"
mid.mid.col = "#fa88d4"
post.post.col = "#6868fc"


#Load full native-space DataFrame
df.orig <- read.table("dataframe.csv")
df.age <- read.table("data-age.csv",
                     sep = ",", header = T) %>%
  separate(Identifiers, into = c("subject","session"), sep = ",") %>%
  dplyr::select(-session)
df.download <- read.table("participants.tsv",
                          sep = "\t", header = T) %>%
  separate(participant_id, into = c("trash","subject"), sep = "-") %>%
  dplyr::select(-trash) %>%
  mutate(sex.d = recode(sex, "M" = 0, "F" = 1, .default = NaN))

df.orig <- merge.data.frame(df.orig, df.download, by = "subject")
df.orig <- merge.data.frame(df.orig, df.age, by = "subject") %>%
  rename(age = "age_04")

#Clean data
df.orig <- df.orig %>%
  filter(na.ratio > 0.7) %>%                #Missing no more than 30% of a given ROI to drop-out
  filter(mean_fd < 0.55) %>%                #mean frame displacement no more than 0.55 mm
  dplyr::group_by(acq,axis,hemi) %>%
  mutate(median_ivs = median(ivs)) %>%      #remove outliers grouped by acquisition, parcel, and hemisphere
  mutate(iqr_ivs = IQR(ivs)) %>%
  filter(!(abs(ivs - median_ivs) > 1.5*iqr_ivs)) %>%
  ungroup()

#recode variables
df.orig$acq <- as.factor(df.orig$acq)
df.orig$hemi <- as.factor(df.orig$hemi)
df.orig$axis <- as.factor(df.orig$axis)
df.orig$age.f <- as.factor(df.orig$age > 35)
df.orig$age.f <- recode(df.orig$age.f, "FALSE" = "young", "TRUE" = "old")
df.orig$age.d <- recode(df.orig$age.f, "young" = 0, "old" = 1)
df.orig$subject <- as.factor(df.orig$subject)

############################################
########    Proportional Parcels    ########
############################################


#############
####1400#####
#############

#Terciles
df.tercile.1400 <- df.orig %>% 
  filter(age.d == 0) %>%
  filter(acq == "1400") %>%
  filter(axis == "ant" | 
           axis == "mid" |
           axis == "mid+post" |
           axis == "post") %>%
  mutate(axis = recode(axis,
                       "ant" = "Anterior",
                       "mid" = "Middle",
                       "mid+post" = "Middle + Posterior",
                       "post" = "Posterior"))

df.dems.1400 <- df.tercile.1400 %>%
  distinct(subject, age, sex)

df.tercile.1400$num_voxels.c <- scale(df.tercile.1400$num_voxels, center = T, scale = F)
df.tercile.1400$mean_dist.c <- scale(df.tercile.1400$mean_dist, center = T, scale = F)
df.tercile.1400$x_dist.c <- scale(df.tercile.1400$x_dist, center = T, scale = F)
df.tercile.1400$y_dist.c <- scale(df.tercile.1400$y_dist, center = T, scale = F)
df.tercile.1400$z_dist.c <- scale(df.tercile.1400$z_dist, center = T, scale = F)
df.tercile.1400$fix_ratio.c <- scale(df.tercile.1400$fix_ratio, center = T, scale = F)
df.tercile.1400$mean_fd.c <- scale(df.tercile.1400$mean_fd, center = T, scale = F)
df.tercile.1400$age.c <- scale(df.tercile.1400$age, center = T, scale = F)
df.tercile.1400$sex.e <- df.tercile.1400$sex.d - 0.5
df.tercile.1400$hemi.e <- df.tercile.1400$hemi.d - 0.5


final_subjs.1400 <- unique(df.tercile.1400$subject)
write.table(final_subjs.1400, "final_subjs_1400.txt", quote = F, row.names = F, col.names = F)

#number of voxels plot
plt.terciles.1400.voxels <- ggplot(data=df.tercile.1400, mapping=aes(acq, num_voxels, color = axis)) +
  geom_violin(position = dodgewhole, alpha = 0, size = 0.5) +
   ylab("voxels") + labs(fill = "parcel", color = "parcel") +
  viridis::scale_color_viridis(discrete = 4, direction = -1)


plt.terciles.1400.voxels

#mean distance plot
plt.terciles.1400.mean_dist <- ggplot(data=df.tercile.1400, mapping=aes(hemi,mean_dist, color = axis)) +
  geom_violin(position = dodgewhole, alpha = 0, size = 0.5) +
  xlab("hemi") + ylab("mean distance from center (mm)") + labs(fill = "parcel", color = "parcel") +
  viridis::scale_color_viridis(discrete = 4, direction = -1)

plt.terciles.1400.mean_dist

#####
## No Covariates
####
mdl.tercile.1400.nocov.f <- lme4::lmer(ivs ~ hemi*axis + (1|subject), df.tercile.1400)
mdl.tercile.1400.nocov.e <- lme4::lmer(ivs ~ hemi.e*axis + (1|subject), df.tercile.1400)
emms.tercile.1400.nocov <- emmeans(mdl.tercile.1400.nocov.e, (pairwise~axis))
p.con.tercile.1400.nocov <- contrast(emms.tercile.1400.nocov[[1]], interaction = c("pairwise"), by = NULL, adjust = "tukey")
#p.concon.tercile.1400.nocov <- contrast(emms.tercile.1400.nocov[[2]],interaction = c("pairwise"), adjust = "tukey")

#plot_hemispheres
lsm.tercile.1400.nocov.f.hemis <- emmeans(mdl.tercile.1400.nocov.e,(pairwise~axis))
lsm.tercile.1400.nocov.hemis <- lsm.tercile.1400.nocov.f.hemis$emmeans[1:4,]

## binary plot
lsm.binary.1400.nocov.hemis <- as.data.frame(lsm.tercile.1400.nocov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

df.binary.1400 <- df.tercile.1400 %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

rownames(lsm.binary.1400.nocov.hemis) <- NULL
lsm.binary.1400.nocov.hemis$axis <- as.character(c('Anterior','Middle + Posterior'))

plt.binary.1400.nocov.hemis <- ggplot(data=lsm.binary.1400.nocov.hemis, 
                                  mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.binary.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.binary.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid_post.col)) +
  scale_fill_manual(values = c(ant.col,mid_post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.binary.1400.nocov.hemis

ggsave('plots/1400-binary-nocov.png',plt.binary.1400.nocov.hemis, width = 5, height = 5)

 ## ternary plot
lsm.ternary.1400.nocov.hemis <- as.data.frame(lsm.tercile.1400.nocov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

df.ternary.1400 <- df.tercile.1400 %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

rownames(lsm.ternary.1400.nocov.hemis) <- NULL
lsm.ternary.1400.nocov.hemis$axis <- as.character(c('Anterior','Middle', 'Posterior'))

plt.ternary.1400.nocov.hemis <- ggplot(data=lsm.ternary.1400.nocov.hemis, 
                                      mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.ternary.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.ternary.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid.col,post.col)) +
  scale_fill_manual(values = c(ant.col,mid.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.ternary.1400.nocov.hemis

ggsave('plots/1400-ternary-nocov.png',plt.ternary.1400.nocov.hemis,width = 5, height = 5)

## combo plot
lsm.combo.1400.nocov.hemis <- as.data.frame(lsm.tercile.1400.nocov.hemis) %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

df.combo.1400 <- df.tercile.1400 %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

rownames(lsm.combo.1400.nocov.hemis) <- NULL
lsm.combo.1400.nocov.hemis$axis <- as.character(c('Middle','Middle + Posterior', 'Posterior'))

plt.combo.1400.nocov.hemis <- ggplot(data=lsm.combo.1400.nocov.hemis, 
                                       mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.combo.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.combo.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(mid.col,mid_post.col,post.col)) +
  scale_fill_manual(values = c(mid.col,mid_post.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.combo.1400.nocov.hemis

ggsave('plots/1400-combo-nocov.png',plt.combo.1400.nocov.hemis,width = 5, height = 5)


#####
## Distance Covariates
####
mdl.tercile.1400.cov.e <- lme4::lmer(ivs ~ hemi.e*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.tercile.1400)

emms.tercile.1400.cov <- emmeans(mdl.tercile.1400.cov.e, (pairwise~axis))
p.con.tercile.1400.cov <- contrast(emms.tercile.1400.cov[[1]], interaction = c("pairwise"), by = NULL, adjust = "tukey")

#plot_hemispheres
lsm.tercile.1400.cov.f.hemis <- emmeans(mdl.tercile.1400.cov.e,(pairwise~axis))
lsm.tercile.1400.cov.hemis <- lsm.tercile.1400.cov.f.hemis$emmeans[1:4,]
lsm.tercile.1400.cov.hemis <- as.data.frame(lsm.tercile.1400.cov.hemis)

rownames(lsm.tercile.1400.cov.hemis) <- NULL
lsm.tercile.1400.cov.hemis$axis <- as.character(c('Anterior','Middle','Middle + Posterior','Posterior'))

## binary plot
lsm.binary.1400.cov.hemis <- as.data.frame(lsm.tercile.1400.cov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

df.binary.1400 <- df.tercile.1400 %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

rownames(lsm.binary.1400.cov.hemis) <- NULL
lsm.binary.1400.cov.hemis$axis <- as.character(c('Anterior','Middle + Posterior'))

plt.binary.1400.cov.hemis <- ggplot(data=lsm.binary.1400.cov.hemis, 
                                      mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.binary.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.binary.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid_post.col)) +
  scale_fill_manual(values = c(ant.col,mid_post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.binary.1400.cov.hemis

ggsave('plots/1400-binary-cov.png',plt.binary.1400.cov.hemis, width = 5, height = 5)

## ternary plot
lsm.ternary.1400.cov.hemis <- as.data.frame(lsm.tercile.1400.cov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

df.ternary.1400 <- df.tercile.1400 %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

rownames(lsm.ternary.1400.cov.hemis) <- NULL
lsm.ternary.1400.cov.hemis$axis <- as.character(c('Anterior','Middle', 'Posterior'))

plt.ternary.1400.cov.hemis <- ggplot(data=lsm.ternary.1400.cov.hemis, 
                                       mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.ternary.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.ternary.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid.col,post.col)) +
  scale_fill_manual(values = c(ant.col,mid.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.ternary.1400.cov.hemis

ggsave('plots/1400-ternary-cov.png',plt.ternary.1400.cov.hemis,width = 5, height = 5)

## combo plot
lsm.combo.1400.cov.hemis <- as.data.frame(lsm.tercile.1400.cov.hemis) %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

df.combo.1400 <- df.tercile.1400 %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

rownames(lsm.combo.1400.cov.hemis) <- NULL
lsm.combo.1400.cov.hemis$axis <- as.character(c('Middle','Middle + Posterior', 'Posterior'))

plt.combo.1400.cov.hemis <- ggplot(data=lsm.combo.1400.cov.hemis, 
                                     mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.combo.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.combo.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(mid.col,mid_post.col,post.col)) +
  scale_fill_manual(values = c(mid.col,mid_post.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(0,0.025)

plt.combo.1400.cov.hemis

ggsave('plots/1400-combo-cov.png',plt.combo.1400.cov.hemis, width = 5, height = 5)


# Sixths
df.sixths.1400 <- df.orig %>% 
  filter(age.d == 0) %>%
  filter(acq == "1400") %>%
  filter(axis == "one" | 
           axis == "two" |
           axis == "three" |
           axis == "four" |
           axis == "five" |
           axis == "six") %>%
  mutate(axis = factor(axis, levels = c('one','two','three','four','five','six'))) %>%
  mutate(axis.e = recode(axis,
                         "one" = 1,
                         "two" = 2,
                         "three" = 3,
                         "four" = 4,
                         "five" = 5,
                         "six" = 6)) %>%
  mutate(axis.c = scale(axis.e, center = T, scale = F))

df.sixths.1400$num_voxels.c <- as.vector(scale(df.sixths.1400$num_voxels, center = T, scale = F))
df.sixths.1400$mean_dist.c <- as.vector(scale(df.sixths.1400$mean_dist, center = T, scale = F))
df.sixths.1400$axis.e.quad <- as.vector(df.sixths.1400$axis.e^2)
df.sixths.1400$x_dist.c <- as.vector(scale(df.sixths.1400$x_dist, center = T, scale = F))
df.sixths.1400$y_dist.c <- as.vector(scale(df.sixths.1400$y_dist, center = T, scale = F))
df.sixths.1400$z_dist.c <- as.vector(scale(df.sixths.1400$z_dist, center = T, scale = F))
df.sixths.1400$mean_fd.c <- as.vector(scale(df.sixths.1400$mean_fd, center = T, scale = F))
df.sixths.1400$age.c <- as.vector(scale(df.sixths.1400$age, center = T, scale = F))
df.sixths.1400$hemi.e <- df.sixths.1400$hemi.d - 0.5
df.sixths.1400$sex.e <- df.sixths.1400$sex.d - 0.5
df.sixths.1400$hemi.rev <- factor(df.sixths.1400$hemi, levels = c("r","l"))

#MODEL HEMI
mdl.sixths.1400.hemis <- lme4::lmer(ivs ~ hemi.e*axis.e + hemi.e*axis.e.quad + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.1400)
summary(mdl.sixths.1400.hemis)


#MODEL for stats
mdl.sixths.1400.quad.forstats <- lmer(ivs ~ hemi.e*poly(axis.c,2) + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.1400)
mdl.sixths.1400.linear.forstats <- lmer(ivs ~ hemi.e*poly(axis.c,1) + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.1400)
summ.sixths.1400 <- anova(mdl.sixths.1400.quad.forstats, mdl.sixths.1400.linear.forstats)

betas.sixths.1400 <- summary(mdl.sixths.1400.quad.forstats, ddf = "Kenward-Roger")

df.sixths.1400.predict.hemis <- data.frame(
  crossing(hemi.e = c(0),
           axis.e = seq(1,6,by = 0.01)),
  x_dist.c = as.vector(0),
  y_dist.c = as.vector(0),
  z_dist.c = as.vector(0)
)

df.sixths.1400.predict.hemis$axis.e.quad <- df.sixths.1400.predict.hemis$axis.e^2

mat <- model.matrix(~ hemi.e*axis.e + hemi.e*axis.e.quad + x_dist.c + y_dist.c + z_dist.c , df.sixths.1400.predict.hemis)
df.sixths.1400.predict.hemis$fit <- predict(mdl.sixths.1400.hemis, df.sixths.1400.predict.hemis, re.form=NA)
predscov <- diag(mat %*% tcrossprod(vcov(mdl.sixths.1400.hemis),mat))
df.sixths.1400.predict.hemis$lwr <- df.sixths.1400.predict.hemis$fit-1.96*sqrt(predscov)
df.sixths.1400.predict.hemis$upr <- df.sixths.1400.predict.hemis$fit+1.96*sqrt(predscov)


## Factor model
mdl.sixths.1400.f <- lme4::lmer(ivs ~ hemi.e*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.1400)

lsm.sixths.1400.f.hemis <- emmeans(mdl.sixths.1400.f,(pairwise~axis))
lsm.sixths.1400.hemis <- lsm.sixths.1400.f.hemis$emmeans[1:6,]
lsm.sixths.1400.hemis <- as.data.frame(lsm.sixths.1400.hemis)

rownames(lsm.sixths.1400.hemis) <- NULL
lsm.sixths.1400.hemis$axis <- factor(c('one','two','three','four','five','six'),
                                     levels = c('one','two','three','four','five','six'))

plt.sixths.1400.hemis <- ggplot(data=lsm.sixths.1400.hemis, 
                                mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.sixths.1400, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.2) +
  geom_point(data=df.sixths.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 0.8), alpha = 0.5, size = 0.4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.24), size = 1, width = 0, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.24), size = 2)+
  ylab("inter-voxel similarity") + 
  labs(color = "parcel", fill = "parcel") +
  #viridis::scale_color_viridis(discrete = 4, direction = -1) +
  scale_color_manual(values = c(ant.col,ant.mid.col,mid.col,mid.mid.col,post.post.col,post.col)) +
  #viridis::scale_fill_viridis(discrete = 4, direction = -1) +
  scale_fill_manual(values = c(ant.col,ant.mid.col,mid.col,mid.mid.col,post.post.col,post.col)) +
  geom_line(data = df.sixths.1400.predict.hemis, mapping=aes(axis.e, fit), color = "black") +
  geom_ribbon(df.sixths.1400.predict.hemis, mapping=aes(axis.e, ymax = upr, ymin = lwr), 
              color = "gray", fill = "gray", alpha = 0.4) +
  theme(axis.text.x=element_blank())


plt.sixths.1400.hemis

ggsave('plots/1400-sixths-cov.png',plt.sixths.1400.hemis, width = 5, height = 8)

#############
####CAP#####
#############

#Terciles
df.tercile.CAP <- df.orig %>% 
  filter(age.d == 0) %>%
  filter(acq == "CAP") %>%
  filter(axis == "ant" | 
           axis == "mid" |
           axis == "mid+post" |
           axis == "post") %>%
  mutate(axis = recode(axis,
                       "ant" = "Anterior",
                       "mid" = "Middle",
                       "mid+post" = "Middle + Posterior",
                       "post" = "Posterior"))

df.dems.CAP <- df.tercile.CAP %>%
  distinct(subject, age, sex)

df.tercile.CAP$num_voxels.c <- scale(df.tercile.CAP$num_voxels, center = T, scale = F)
df.tercile.CAP$mean_dist.c <- scale(df.tercile.CAP$mean_dist, center = T, scale = F)
df.tercile.CAP$x_dist.c <- scale(df.tercile.CAP$x_dist, center = T, scale = F)
df.tercile.CAP$y_dist.c <- scale(df.tercile.CAP$y_dist, center = T, scale = F)
df.tercile.CAP$z_dist.c <- scale(df.tercile.CAP$z_dist, center = T, scale = F)
df.tercile.CAP$fix_ratio.c <- scale(df.tercile.CAP$fix_ratio, center = T, scale = F)
df.tercile.CAP$mean_fd.c <- scale(df.tercile.CAP$mean_fd, center = T, scale = F)
df.tercile.CAP$age.c <- scale(df.tercile.CAP$age, center = T, scale = F)
df.tercile.CAP$sex.e <- df.tercile.CAP$sex.d - 0.5
df.tercile.CAP$hemi.e <- df.tercile.CAP$hemi.d - 0.5


final_subjs.CAP <- unique(df.tercile.CAP$subject)
write.table(final_subjs.CAP, "final_subjs_CAP.txt", quote = F, row.names = F, col.names = F)

#number of voxels plot
plt.terciles.CAP.voxels <- ggplot(data=df.tercile.CAP, mapping=aes(acq, num_voxels, color = axis)) +
  geom_violin(position = dodgewhole, alpha = 0, size = 0.5) +
   ylab("voxels") + labs(fill = "parcel", color = "parcel") +
  viridis::scale_color_viridis(discrete = 4, direction = -1)


plt.terciles.CAP.voxels

#mean distance plot
plt.terciles.CAP.mean_dist <- ggplot(data=df.tercile.CAP, mapping=aes(hemi,mean_dist, color = axis)) +
  geom_violin(position = dodgewhole, alpha = 0, size = 0.5) +
  xlab("hemi") + ylab("mean distance from center (mm)") + labs(fill = "parcel", color = "parcel") +
  viridis::scale_color_viridis(discrete = 4, direction = -1)

plt.terciles.CAP.mean_dist

#####
## No Covariates
####
mdl.tercile.CAP.nocov.f <- lme4::lmer(ivs ~ hemi*axis + (1|subject), df.tercile.CAP)
mdl.tercile.CAP.nocov.e <- lme4::lmer(ivs ~ hemi.e*axis + (1|subject), df.tercile.CAP)
emms.tercile.CAP.nocov <- emmeans(mdl.tercile.CAP.nocov.e, (pairwise~axis))
p.con.tercile.CAP.nocov <- contrast(emms.tercile.CAP.nocov[[1]], interaction = c("pairwise"), by = NULL, adjust = "tukey")
#p.concon.tercile.CAP.nocov <- contrast(emms.tercile.CAP.nocov[[2]],interaction = c("pairwise"), adjust = "tukey")

#plot_hemispheres
lsm.tercile.CAP.nocov.f.hemis <- emmeans(mdl.tercile.CAP.nocov.e,(pairwise~axis))
lsm.tercile.CAP.nocov.hemis <- lsm.tercile.CAP.nocov.f.hemis$emmeans[1:4,]

## binary plot
lsm.binary.CAP.nocov.hemis <- as.data.frame(lsm.tercile.CAP.nocov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

df.binary.CAP <- df.tercile.CAP %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

rownames(lsm.binary.CAP.nocov.hemis) <- NULL
lsm.binary.CAP.nocov.hemis$axis <- as.character(c('Anterior','Middle + Posterior'))

plt.binary.CAP.nocov.hemis <- ggplot(data=lsm.binary.CAP.nocov.hemis, 
                                      mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.binary.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.binary.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid_post.col)) +
  scale_fill_manual(values = c(ant.col,mid_post.col)) +
  theme(axis.text.x=element_blank()) + ylim(-0.01, 0.07)

plt.binary.CAP.nocov.hemis

ggsave('plots/CAP-binary-nocov.png',plt.binary.CAP.nocov.hemis, width = 5, height = 5)

## ternary plot
lsm.ternary.CAP.nocov.hemis <- as.data.frame(lsm.tercile.CAP.nocov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

df.ternary.CAP <- df.tercile.CAP %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

rownames(lsm.ternary.CAP.nocov.hemis) <- NULL
lsm.ternary.CAP.nocov.hemis$axis <- as.character(c('Anterior','Middle', 'Posterior'))

plt.ternary.CAP.nocov.hemis <- ggplot(data=lsm.ternary.CAP.nocov.hemis, 
                                       mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.ternary.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.ternary.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid.col,post.col)) +
  scale_fill_manual(values = c(ant.col,mid.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(-0.01,0.07)

plt.ternary.CAP.nocov.hemis

ggsave('plots/CAP-ternary-nocov.png',plt.ternary.CAP.nocov.hemis, width = 5, height = 5)

## combo plot
lsm.combo.CAP.nocov.hemis <- as.data.frame(lsm.tercile.CAP.nocov.hemis) %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

df.combo.CAP <- df.tercile.CAP %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

rownames(lsm.combo.CAP.nocov.hemis) <- NULL
lsm.combo.CAP.nocov.hemis$axis <- as.character(c('Middle','Middle + Posterior', 'Posterior'))

plt.combo.CAP.nocov.hemis <- ggplot(data=lsm.combo.CAP.nocov.hemis, 
                                     mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.combo.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.combo.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(mid.col,mid_post.col,post.col)) +
  scale_fill_manual(values = c(mid.col,mid_post.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(-0.01,0.07)

plt.combo.CAP.nocov.hemis

ggsave('plots/CAP-combo-nocov.png',plt.combo.CAP.nocov.hemis, width = 5, height = 5)


#####
## Distance Covariates
####
mdl.tercile.CAP.cov.e <- lme4::lmer(ivs ~ hemi.e*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.tercile.CAP)

emms.tercile.CAP.cov <- emmeans(mdl.tercile.CAP.cov.e, (pairwise~axis))
p.con.tercile.CAP.cov <- contrast(emms.tercile.CAP.cov[[1]], interaction = c("pairwise"), by = NULL, adjust = "tukey")

#plot_hemispheres
lsm.tercile.CAP.cov.f.hemis <- emmeans(mdl.tercile.CAP.cov.e,(pairwise~axis))
lsm.tercile.CAP.cov.hemis <- lsm.tercile.CAP.cov.f.hemis$emmeans[1:4,]
lsm.tercile.CAP.cov.hemis <- as.data.frame(lsm.tercile.CAP.cov.hemis)

rownames(lsm.tercile.CAP.cov.hemis) <- NULL
lsm.tercile.CAP.cov.hemis$axis <- as.character(c('Anterior','Middle','Middle + Posterior','Posterior'))

## binary plot
lsm.binary.CAP.cov.hemis <- as.data.frame(lsm.tercile.CAP.cov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

df.binary.CAP <- df.tercile.CAP %>%
  filter(axis == "Anterior" | axis == "Middle + Posterior")

rownames(lsm.binary.CAP.cov.hemis) <- NULL
lsm.binary.CAP.cov.hemis$axis <- as.character(c('Anterior','Middle + Posterior'))

plt.binary.CAP.cov.hemis <- ggplot(data=lsm.binary.CAP.cov.hemis, 
                                    mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.binary.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.binary.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid_post.col)) +
  scale_fill_manual(values = c(ant.col,mid_post.col)) +
  theme(axis.text.x=element_blank())+ ylim(-0.01,0.07)

plt.binary.CAP.cov.hemis

ggsave('plots/CAP-binary-cov.png',plt.binary.CAP.cov.hemis, width = 5, height = 5)

## ternary plot
lsm.ternary.CAP.cov.hemis <- as.data.frame(lsm.tercile.CAP.cov.hemis) %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

df.ternary.CAP <- df.tercile.CAP %>%
  filter(axis == "Anterior" | axis == "Middle" | axis == "Posterior")

rownames(lsm.ternary.CAP.cov.hemis) <- NULL
lsm.ternary.CAP.cov.hemis$axis <- as.character(c('Anterior','Middle', 'Posterior'))

plt.ternary.CAP.cov.hemis <- ggplot(data=lsm.ternary.CAP.cov.hemis, 
                                     mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.ternary.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.ternary.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(ant.col,mid.col,post.col)) +
  scale_fill_manual(values = c(ant.col,mid.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(-0.01,0.07)

plt.ternary.CAP.cov.hemis

ggsave('plots/CAP-ternary-cov.png',plt.ternary.CAP.cov.hemis, width = 5, height = 5)

## combo plot
lsm.combo.CAP.cov.hemis <- as.data.frame(lsm.tercile.CAP.cov.hemis) %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

df.combo.CAP <- df.tercile.CAP %>%
  filter(axis == "Middle" | axis == "Middle + Posterior" | axis == "Posterior")

rownames(lsm.combo.CAP.cov.hemis) <- NULL
lsm.combo.CAP.cov.hemis$axis <- as.character(c('Middle','Middle + Posterior', 'Posterior'))

plt.combo.CAP.cov.hemis <- ggplot(data=lsm.combo.CAP.cov.hemis, 
                                   mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.combo.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.25) +
  geom_point(data=df.combo.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 1), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.15), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.15), size = 3)+
  ylab("inter-voxel similarity") + labs(fill = "parcel", color = "parcel") +
  scale_color_manual(values = c(mid.col,mid_post.col,post.col)) +
  scale_fill_manual(values = c(mid.col,mid_post.col,post.col)) +
  theme(axis.text.x=element_blank()) + ylim(-0.01,0.07)

plt.combo.CAP.cov.hemis

ggsave('plots/CAP-combo-cov.png',plt.combo.CAP.cov.hemis, width = 5, height = 5)


# Sixths
df.sixths.CAP <- df.orig %>% 
  filter(age.d == 0) %>%
  filter(acq == "CAP") %>%
  filter(axis == "one" | 
           axis == "two" |
           axis == "three" |
           axis == "four" |
           axis == "five" |
           axis == "six") %>%
  mutate(axis = factor(axis, levels = c('one','two','three','four','five','six'))) %>%
  mutate(axis.e = recode(axis,
                         "one" = 1,
                         "two" = 2,
                         "three" = 3,
                         "four" = 4,
                         "five" = 5,
                         "six" = 6)) %>%
  mutate(axis.c = scale(axis.e, center = T, scale = F))

df.sixths.CAP$axis.e.quad <- as.vector(df.sixths.CAP$axis.e^2)
df.sixths.CAP$num_voxels.c <- as.vector(scale(df.sixths.CAP$num_voxels, center = T, scale = F))
df.sixths.CAP$mean_dist.c <- as.vector(scale(df.sixths.CAP$mean_dist, center = T, scale = F))
df.sixths.CAP$x_dist.c <- as.vector(scale(df.sixths.CAP$x_dist, center = T, scale = F))
df.sixths.CAP$y_dist.c <- as.vector(scale(df.sixths.CAP$y_dist, center = T, scale = F))
df.sixths.CAP$z_dist.c <- as.vector(scale(df.sixths.CAP$z_dist, center = T, scale = F))
df.sixths.CAP$mean_fd.c <- as.vector(scale(df.sixths.CAP$mean_fd, center = T, scale = F))
df.sixths.CAP$age.c <- as.vector(scale(df.sixths.CAP$age, center = T, scale = F))
df.sixths.CAP$hemi.e <- df.sixths.CAP$hemi.d - 0.5
df.sixths.CAP$sex.e <- df.sixths.CAP$sex.d - 0.5
df.sixths.CAP$hemi.rev <- factor(df.sixths.CAP$hemi, levels = c("r","l"))

#MODEL HEMI
mdl.sixths.CAP.hemis <- lme4::lmer(ivs ~ hemi.e*axis.e + hemi.e*axis.e.quad+ x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.CAP)
summary(mdl.sixths.CAP.hemis)


#MODEL for stats
mdl.sixths.CAP.quad.forstats <- lmer(ivs ~ hemi.e*poly(axis.c,2) + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.CAP)
mdl.sixths.CAP.linear.forstats <- lmer(ivs ~ hemi.e*axis.c + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.CAP)
summ.sixths.CAP <- anova(mdl.sixths.CAP.quad.forstats, mdl.sixths.CAP.linear.forstats)
betas.sixths.CAP <- summary(mdl.sixths.CAP.quad.forstats, ddf = "Kenward-Roger")


df.sixths.CAP.predict.hemis <- data.frame(
  crossing(hemi.e = c(0),
           axis.e = seq(1,6,by = 0.01)),
  x_dist.c = as.vector(0),
  y_dist.c = as.vector(0),
  z_dist.c = as.vector(0)
)
df.sixths.CAP.predict.hemis$axis.e.quad <- df.sixths.CAP.predict.hemis$axis.e^2


mat <- model.matrix(~ hemi.e*axis.e + hemi.e*axis.e.quad + x_dist.c + y_dist.c + z_dist.c , df.sixths.CAP.predict.hemis)
df.sixths.CAP.predict.hemis$fit <- predict(mdl.sixths.CAP.hemis, df.sixths.CAP.predict.hemis, re.form=NA)
predscov <- diag(mat %*% tcrossprod(vcov(mdl.sixths.CAP.hemis),mat))
df.sixths.CAP.predict.hemis$lwr <- df.sixths.CAP.predict.hemis$fit-1.96*sqrt(predscov)
df.sixths.CAP.predict.hemis$upr <- df.sixths.CAP.predict.hemis$fit+1.96*sqrt(predscov)


## Factor model
mdl.sixths.CAP.f <- lme4::lmer(ivs ~ hemi.e*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), df.sixths.CAP)

lsm.sixths.CAP.f.hemis <- emmeans(mdl.sixths.CAP.f,(pairwise~axis))
lsm.sixths.CAP.hemis <- lsm.sixths.CAP.f.hemis$emmeans[1:6,]
lsm.sixths.CAP.hemis <- as.data.frame(lsm.sixths.CAP.hemis)

rownames(lsm.sixths.CAP.hemis) <- NULL
lsm.sixths.CAP.hemis$axis <- factor(c('one','two','three','four','five','six'),
                                     levels = c('one','two','three','four','five','six'))

plt.sixths.CAP.hemis <- ggplot(data=lsm.sixths.CAP.hemis, 
                                mapping=aes(axis, color = axis)) +
  geom_flat_violin(data=df.sixths.CAP, 
                   mapping=aes(axis, ivs, fill = axis), color = "gray", 
                   position = position_nudge(0.1), alpha = 0.5, size = 0.2) +
  geom_point(data=df.sixths.CAP, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 0.8), alpha = 0.5, size = 0.4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.24), size = 1, width = 0, color = "black")+
  geom_point(aes(axis, y = emmean), stat = "identity", position = position_nudge(-0.24), size = 2)+
  ylab("inter-voxel similarity") + 
  labs(color = "parcel", fill = "parcel") +
  #viridis::scale_color_viridis(discrete = 4, direction = -1) +
  scale_color_manual(values = c(ant.col,ant.mid.col,mid.col,mid.mid.col,post.post.col,post.col)) +
  #viridis::scale_fill_viridis(discrete = 4, direction = -1) +
  scale_fill_manual(values = c(ant.col,ant.mid.col,mid.col,mid.mid.col,post.post.col,post.col)) +
  geom_line(data = df.sixths.CAP.predict.hemis, mapping=aes(axis.e, fit), color = "black") +
  geom_ribbon(df.sixths.CAP.predict.hemis, mapping=aes(axis.e, ymax = upr, ymin = lwr), 
              color = "gray", fill = "gray", alpha = 0.4) +
  theme(axis.text.x=element_blank())


plt.sixths.CAP.hemis

ggsave('plots/CAP-sixths-cov.png',plt.sixths.CAP.hemis, width = 5, height = 8)



############################################
#############    Components   ##############
############################################


df.native <- read.table("dataframe.csv")
df.age <- read.table("data-age.csv",
                     sep = ",", header = T) %>%
  separate(Identifiers, into = c("subject","session"), sep = ",") %>%
  dplyr::select(-session)
df.download <- read.table("participants.tsv",
                          sep = "\t", header = T) %>%
  separate(participant_id, into = c("trash","subject"), sep = "-") %>%
  dplyr::select(-trash) %>%
  mutate(sex.d = recode(sex, "M" = 0, "F" = 1, .default = NaN))

df.native <- merge.data.frame(df.native, df.download, by = "subject")
df.native <- merge.data.frame(df.native, df.age, by = "subject") %>%
  rename(age = "age_04")
df.native$acq <- as.factor(df.native$acq)
df.native$hemi <- as.factor(df.native$hemi)
df.native$axis <- as.factor(df.native$axis)

df.native$age.f <- as.factor(df.native$age > 35)
df.native$age.f <- recode(df.native$age.f, "FALSE" = "young", "TRUE" = "old")
df.native$age.d <- recode(df.native$age.f, "young" = 0, "old" = 1)
df.native$axis <- recode(df.native$axis, "mid+post" = "mid.post")
df.native$subject <- as.factor(df.native$subject)

df.native.filt <- df.native %>%
  filter(na.ratio > 0.7) %>%
  filter(mean_fd < 0.55) %>%
  filter(overlap_ratio > 0.7) %>%                     # filter out components that don't overlap the freesurfer mask by more than 20%
  dplyr::group_by(acq,axis,hemi) %>%                  # remove outliers grouped by acq, parcel, and hemisphere
  mutate(median_ivs = median(ivs)) %>%
  mutate(iqr_ivs = IQR(ivs)) %>%
  filter(!(abs(ivs - median_ivs) > 1.5*iqr_ivs)) %>%
  ungroup()



######### Main Model ############

df.blessing.comps.1400 <- df.native.filt %>% 
  filter(age.d == 0) %>%
  filter(acq == "1400") %>%
  filter(axis == "a" |
           axis == "am" |
           axis == "al" |
           axis == "m" |
           axis == "p") %>%
  mutate(axis = factor(axis, levels = c("am","a","al","m","p"))) %>%
  mutate(axis = recode(axis,
                       "a" = "AL",
                       "am" = "AM",
                       "al" = "PAL",
                       "m" = "M",
                       "p" = "P")) %>%
  mutate(hemi = factor(hemi, levels = c("r","l")))

df.dems.blessing.1400 <- df.blessing.comps.1400 %>%
  distinct(subject, age, sex)

df.blessing.comps.1400$num_voxels.c <- scale(df.blessing.comps.1400$num_voxels, center = T, scale = F)
df.blessing.comps.1400$mean_dist.c <- scale(df.blessing.comps.1400$mean_dist, center = T, scale = F)
df.blessing.comps.1400$x_dist.c <- scale(df.blessing.comps.1400$x_dist, center = T, scale = F)
df.blessing.comps.1400$y_dist.c <- scale(df.blessing.comps.1400$y_dist, center = T, scale = F)
df.blessing.comps.1400$z_dist.c <- scale(df.blessing.comps.1400$z_dist, center = T, scale = F)
df.blessing.comps.1400$mean_fd.c <- scale(df.blessing.comps.1400$mean_fd, center = T, scale = F)
df.blessing.comps.1400$age.c <- scale(df.blessing.comps.1400$age, center = T, scale = F)
df.blessing.comps.1400$overlap_ratio.c <- scale(df.blessing.comps.1400$overlap_ratio, center = T, scale = F)
df.blessing.comps.1400$hemi.e <- df.blessing.comps.1400$hemi.d - 0.5


#factor model
mdl.blessing.comps.1400.f <- lme4::lmer(ivs ~ hemi*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), 
                                        df.blessing.comps.1400)
mdl.blessing.comps.1400.e <- lme4::lmer(ivs ~ hemi.e*axis + x_dist.c + y_dist.c + z_dist.c + (1|subject), 
                                        df.blessing.comps.1400)
emms.blessing.1400 <- emmeans(mdl.blessing.comps.1400.f, (pairwise~axis), by = "hemi")
p.con.blessing.1400 <- contrast(emms.blessing.1400[[1]], interaction = c("pairwise"), by = NULL, adjust = "tukey")

## Hemispheres Plot
lsm.blessing.comps.1400.f.hemis <- emmeans(mdl.blessing.comps.1400.f,(pairwise~axis), by = "hemi")
dodgehalf <- position_dodge(width=0.5)
dodgewhole <- position_dodge(width = 1)
lsm.blessing.comps.1400.hemis <- lsm.blessing.comps.1400.f.hemis$emmeans[1:10,]
lsm.blessing.comps.1400.hemis <- as.data.frame(lsm.blessing.comps.1400.hemis)

rownames(lsm.blessing.comps.1400.hemis) <- NULL
lsm.blessing.comps.1400.hemis$axis <- factor(rep(c('AM','AL','PAL',
                                                   'M','P'),times = 2),
                                             levels = c('AM','AL','PAL',
                                                        'M','P'))
lsm.blessing.comps.1400.hemis$hemi <- as.factor(rep(c('l','r'),each = 5))

plt.blessing.comps.1400.hemis <- ggplot(data=lsm.blessing.comps.1400.hemis, 
                                        mapping=aes(axis, fill = axis)) +
  facet_grid(cols=rev(vars(hemi)),labeller = hemisphere_labeller) +
  geom_flat_violin(data=df.blessing.comps.1400, 
                   mapping=aes(axis, ivs, fill = axis), 
                   position = position_nudge(0.1), alpha = 0.7, size = 0.25, color = "gray") +
  geom_point(data=df.blessing.comps.1400, mapping = aes(axis, ivs, color=axis), 
             position = position_jitterdodge(jitter.width = 0.5), alpha = 0.5, size = 0.4) +
  #geom_bar(aes(axis, y = emmean, fill = axis), stat = "identity", position = position_nudge(-0.24), width = 0.28, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_nudge(-0.20), size = 1, width = 0.00, color = "black")+
  geom_point(aes(axis, y = emmean, color = axis), stat = "identity", position = position_nudge(-0.20), size = 3)+
  ylab("inter-voxel similarity") + 
  labs(color = "parcel", fill = "parcel") +
  scale_color_manual(values = c("#D32D1F","#FFEA5D","#187321","#68E2D8","#0000BC")) +
  scale_fill_manual(values = c("#D32D1F","#FFEA5D","#187321","#68E2D8","#0000BC")) +
  theme(axis.text.x = element_text(size = 24), axis.title.y = element_text(size = 24))


plt.blessing.comps.1400.hemis

ggsave('plots/1400-blessing-cov.png',plt.blessing.comps.1400.hemis, width = 10, height = 6)
