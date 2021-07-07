#load packages
library(plyr)
library(vip)
library(ranger)
library(effects)
library(ggeffects)
library(MASS)
library(psych)
library(lsr)
library(dplyr)
library(regclass)
library(reshape2)
library(ggplot2)
library(rockchalk)
library(corrplot)
library(lme4)
library(lmerTest)
set.seed(1)


##### FUNCTIONS ####
#function to transform all variables into normal shape 
box_cox_transf <- function(x, min, max) {
  bc = boxcox(x ~ 1, lambda = seq(min, max, 0.1), plotit = FALSE)
  df = data.frame(bc$x, bc$y)
  df2 = df[with(df, order(-df$bc.y)),]
  lambda = df2[1, "bc.x"]
  print(lambda)
  if (lambda != 0) {
    x.transf = (x ^ lambda - 1)/lambda
  } else {
    x.transf = log(x)
  }
  return(x.transf)
}

#function to get upper triangle from correlation matrix 
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

##### DATA PRE-PROCESSING #####
#load the data
data <- read.csv(file = "output/fsc_measures.csv")
iconicity <- read.csv(file = '/media/gioca90/University/TiU/Research/Resources/iconicity_ratings.csv')
df.mald = read.csv(file = '/media/gioca90/University/TiU/Research/Resources/MALD/AllData.txt', sep = '\t')
zeno_aoa <- read.csv(file = "data/Zeno et al Data.csv") 

#transform morph into factor
data$morph <- factor(data$morph,    
                     level=c(0,1),    
                     labels=c("Mono","Poly"))

# melt Zeno et al df so grade becomes a single variable
zeno_aoa <- melt(zeno_aoa, 
                 id.vars = c('Word', 'Iconicity', 'Length', 'Concreteness', 'LgChildDirectedFreq', 'POS', 'OLD'), 
                 measure.vars = c('Gr1', 'Gr2', 'Gr3', 'Gr4', 'Gr5', 'Gr6', 'Gr7', 'Gr8', 'Gr9', 'Gr10', 'Gr11', 'Gr12', 'Gr13.'),
                 variable.name = 'Grade', value.name = 'ZenoFrequency')
zeno_aoa$Grade <- gsub('Gr', '', zeno_aoa$Grade)
zeno_aoa$Grade <- as.numeric(zeno_aoa$Grade)

zeno_aoa <- na.omit(zeno_aoa)

zeno_aoa$ZenoFrequency <- log(zeno_aoa$ZenoFrequency + 1)
zeno_aoa$Iconicity.z <- c(scale(box_cox_transf(zeno_aoa$Iconicity + abs(min(zeno_aoa$Iconicity)) + 0.01, -10, 10), center = T, scale = T))


### Transform variables into normal shape and comparable units (SD) ###
df.ortho <- data[,FALSE]
df.phono <- data[,FALSE]

df.ortho$word = data$word
df.phono$word = data$word

df.ortho$morph = data$morph
df.phono$morph = data$morph

df.ortho$concr.z = c(scale(box_cox_transf(data$concr, -10, 10), center = T, scale = T))
df.phono$concr.z = df.ortho$concr.z

df.ortho$val.z = c(scale(box_cox_transf(data$val, -10, 10), center = T, scale = T))
df.phono$val.z = df.ortho$val.z

df.ortho$freq.z = c(scale(box_cox_transf(data$freq, -10, 10), center = T, scale = T))
df.phono$freq.z = df.ortho$freq.z

df.ortho$snd.z = c(scale(box_cox_transf(data$snd, -10, 10), center = T, scale = T))
df.phono$snd.z = df.ortho$snd.z

df.ortho$n_chars.z = c(scale(box_cox_transf(data$length, -10, 10), center = T, scale = T))
df.ortho$old.z = c(scale(box_cox_transf(data$old20, -10, 10), center = T, scale = T))
df.ortho$OSC_te.z = c(scale(box_cox_transf(data$OSC_te + 0.0001, -10, 10), center = T, scale = T))
df.ortho$OSC_ld.z = c(scale(box_cox_transf(data$OSC_ld, -10, 10), center = T, scale = T))

df.ortho = merge(df.ortho, zeno_aoa[c('Word', 'Iconicity.z', 'Grade', 'ZenoFrequency')], by.x='word', by.y='Word')


df.phono$n_phons.z = c(scale(box_cox_transf(data$n_phon, -10, 10), center = T, scale = T))
df.phono$PSC_te.z = c(scale(box_cox_transf(data$PSC_te + 0.0001, -10, 10), center = T, scale = T))
df.phono$PSC_ld.z = c(scale(box_cox_transf(data$PSC_ld, -10, 10), center = T, scale = T))
df.phono = merge(df.phono, df.mald[c('Item', 'PhonND')], by.x='word', by.y='Item')
df.phono = distinct(df.phono)
df.phono$pnd.z = c(scale(box_cox_transf(df.phono$PhonND + 0.0001, -10, 10), center = T, scale = T))

df.phono = merge(df.phono, zeno_aoa[c('Word', 'Iconicity.z', 'Grade', 'ZenoFrequency')], by.x='word', by.y='Word')


##### MULTIPLE REGRESSION - ORTHOGRAPHY #####
# Base Model
base_model.ortho <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*Grade + snd.z*Grade + morph*Grade + Iconicity.z*Grade + (1|word), 
  data = df.ortho
)
# all interactions are significant
summary(base_model.ortho)
base.ortho.aic <- AIC(base_model.ortho)
# AIC: 39368.84

quantile(df.ortho$snd.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
quantile(df.ortho$old.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
quantile(df.ortho$Iconicity.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

base_model.ortho.old.df = as.data.frame(
  ggpredict(
    base_model.ortho, terms = c(
      "old.z[-1.872, -1.622, -1.185, -0.81, -0.641, -0.406, -0.260, -0.122, 0.581]",
      "Grade[1,5,9,13]"
      )
    )
)

base_model.ortho.snd.df = as.data.frame(
  ggpredict(
    base_model.ortho, terms = c(
      "snd.z[-1.740, -1.325, -1.014, -0.75, -0.475, -0.210, 0.054, 0.380, 0.823]",
      "Grade[1,5,9,13]"
    )
  )
)

base_model.ortho.icon.df = as.data.frame(
  ggpredict(
    base_model.ortho, terms = c(
      "Iconicity.z[1.163, -0.831, -0.553, -0.318, -0.081, 0.154, 0.428, 0.702, 1.298]",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(base_model.ortho.icon.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab('Iconicity') +
  ylab('log(Freq)') +
  ggtitle(expression("2-way interaction of Iconicity and Grade on Frequency (ortho)")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


base_model.ortho.snd_old.df = as.data.frame(
  ggpredict(
    base_model.ortho, terms = c(
      "old.z[-1.872, -1.622, -1.185, -0.81, -0.641, -0.406, -0.260, -0.122, 0.581]",
      "snd.z[-2.103, -0.475, 1.128]",
      "Grade[1,5,9,13]")
  )
)
base_model.ortho.snd_old.df$group <- mapvalues(
  base_model.ortho.snd_old.df$group, 
  from = c("-2.103", "-0.475", "1.128"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.ortho.snd_old.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('OLD20') +
  ylab('Freq') +
  ggtitle(expression("Additive effect of OLD20 and SND on Frequency per Grade")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Base Model with snd*old interaction
base_model.ortho.int <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + (1|word), data = df.ortho
)
base_model.ortho.int.aic = AIC(base_model.ortho.int)
# AIC: 39369.45
summary(base_model.ortho.int)
# target 3-way interaction non significant, all other interactions are
base.ortho.aic - base_model.ortho.int.aic
# deltaAIC: -0.6003026

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
quantile(df.ortho$old.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
base_model.ortho.int.df = as.data.frame(
  ggpredict(
    base_model.ortho.int, terms = c(
      "old.z[-1.872, -1.622, -1.185, -0.81, -0.641, -0.406, -0.260, -0.122, 0.581]",
      "snd.z[-2.103, -0.475, 1.128]",
      "Grade[1,5,9,13]")
  )
)
base_model.ortho.int.df$group <- mapvalues(
  base_model.ortho.int.df$group, 
  from = c("-2.103", "-0.475", "1.128"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.ortho.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('OLD20') +
  ylab('Freq') +
  ggtitle(expression("3-way interaction of OLD20, SND and Grade on Frequency per Grade")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# target embedded neighbors
osc_model.te <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_te.z*Grade + (1|word), 
  data = df.ortho
)
osc_model.te.aic = AIC(osc_model.te)
# 39321.86
summary(osc_model.te)
# Grade*OSCte interaction significant
base.ortho.aic - osc_model.te.aic
# 46.98543


quantile(df.ortho$OSC_te.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
osc_model.te.df = as.data.frame(
  ggpredict(
    osc_model.te, terms = c(
      "OSC_te.z[-1.0382228, -0.6417030, -0.1904866, 0.1565135, 0.4072468, 0.6208112, 0.8579672, 1.0633825, 1.3010812]",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(osc_model.te.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('OSC'[te])) +
  ylab('log(Freq)') +
  ggtitle(expression("2-way interaction of OSC[te] and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# OSCte*morph interaction
osc_by_morph_model.te <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_te.z*Grade*morph + (1|word), 
  data = df.ortho
)
osc_by_morph_model.te.aic = AIC(osc_by_morph_model.te)
# 39336.31
summary(osc_by_morph_model.te)
# target interaction non significant
osc_model.te.aic - osc_by_morph_model.te.aic
# -14.45229



# Levenshtein neighbors
osc_model.ld <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_ld.z*Grade + (1|word), 
  data = df.ortho
)
osc_model.ld.aic = AIC(osc_model.ld)
# 38343.15
summary(osc_model.ld)
# target interaction significant
base_model.ortho.int.aic - osc_model.ld.aic
# 1026.297

quantile(df.ortho$OSC_ld.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
osc_model.ld.df = as.data.frame(
  ggpredict(
    osc_model.ld, terms = c(
      "OSC_ld.z[-0.5849087, -0.1101741, 0.1801345, 0.4086641, 0.6000482, 0.7776773, 1.0016594, 1.3132226, 1.8524416]",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(osc_model.ld.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("2-way interaction of OSC[ld] and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# OSCld*morph interaction
osc_by_morph_model.ld <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + Iconicity.z*Grade + OSC_ld.z*Grade*morph + (1|word), 
  data = df.ortho
)
osc_by_morph_model.ld.aic = AIC(osc_by_morph_model.ld)
# 38351.58
summary(osc_by_morph_model.ld)
# target interaction significant, but AIC is worse
osc_model.ld.aic - osc_by_morph_model.ld.aic
# -8.430733

osc_by_morph_model.ld.df = as.data.frame(
  ggpredict(
    osc_by_morph_model.ld, terms = c(
      "OSC_ld.z[-0.5849087, -0.1101741, 0.1801345, 0.4086641, 0.6000482, 0.7776773, 1.0016594, 1.3132226, 1.8524416]",
      "morph",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(osc_by_morph_model.ld.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("3-way interaction of OSC[ld], morphological complexity and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 



##### MULTIPLE REGRESSION - PHONOLOGY #####
# Base Model
base_model.phono <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*Grade + snd.z*Grade + morph*Grade + Iconicity.z*Grade + (1|word), data = df.phono
)
summary(base_model.phono)
# all interactions are significant
base.phono.aic <- AIC(base_model.phono)
# 36737.08

quantile(df.phono$snd.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
quantile(df.phono$pnd.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
quantile(df.phono$Iconicity.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

base_model.phono.pnd.df = as.data.frame(
  ggpredict(
    base_model.phono, terms = c(
      "pnd.z[-0.43899943, -0.08205788, 0.32963842, 0.55784654, 0.75051840, 0.93544343, 1.12332210, 1.35551414, 1.58609120]",
      "Grade[1,5,9,13]"
    )
  )
)

base_model.phono.snd.df = as.data.frame(
  ggpredict(
    base_model.phono, terms = c(
      "snd.z[-1.7428075, -1.3377134, -1.0402352, -0.7704514, -0.4893360, -0.2505143, 0.0327292, 0.3332677, 0.7730227]",
      "Grade[1,5,9,13]"
    )
  )
)

base_model.phono.icon.df = as.data.frame(
  ggpredict(
    base_model.phono, terms = c(
      "Iconicity.z[-1.1629804, -0.8305493, -0.5525022, -0.3182077, -0.1003609, 0.1537658, 0.4107670, 0.7018816, 1.2284187]",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(base_model.phono.pnd.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab('PND') +
  ylab('log(Freq)') +
  ggtitle(expression("2-way interaction of PND and Grade on Frequency (phono)")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


base_model.phono.pnd_snd.df = as.data.frame(
  ggpredict(
    base_model.phono, terms = c(
      "pnd.z[-0.43899943, -0.08205788, 0.32963842, 0.55784654, 0.75051840, 0.93544343, 1.12332210, 1.35551414, 1.58609120]",
      "snd.z[-2.102646, -0.489336, 1.088724]",
      "Grade[1,5,9,13]")
  )
)
base_model.phono.pnd_snd.df$group <- mapvalues(
  base_model.phono.pnd_snd.df$group, 
  from = c("-2.102646", "-0.489336", "1.088724"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.phono.pnd_snd.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('PND') +
  ylab('Freq') +
  ggtitle(expression("Additive effect of PND and SND on Frequency per Grade")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Base Model with snd*pnd interaction
base_model.phono.int <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + (1|word), data = df.phono
)
base_model.phono.int.aic = AIC(base_model.phono.int)
# 36736.75
summary(base_model.phono.int)
# target interaction is significant, but AIC isn't better
base.phono.aic - base_model.phono.int.aic
# 0.3235339

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
quantile(df.phono$pnd.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
base_model.phono.int.df = as.data.frame(
  ggpredict(
    base_model.phono.int, terms = c(
      "pnd.z[-0.43899943, -0.08205788, 0.32963842, 0.55784654, 0.75051840, 0.93544343, 1.12332210, 1.35551414, 1.58609120]",
      "snd.z[-2.102646, -0.489336, 1.088724]",
      "Grade[1,5,9,13]")
  )
)
base_model.phono.int.df$group <- mapvalues(
  base_model.phono.int.df$group, 
  from = c("-2.102646", "-0.489336", "1.088724"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.phono.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('PND') +
  ylab('Freq') +
  ggtitle(expression("3-way interaction of PND and SND on Frequency per Grade")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# target embedded neighbors
PSC_model.te <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + PSC_te.z*Grade + (1|word), 
  data = df.phono
)
PSC_model.te.aic = AIC(PSC_model.te)
# 36752.03
summary(PSC_model.te)
# target interaction isn't significant
base.phono.aic - PSC_model.te.aic
# -14.95783


# PSCte*morph interaction
PSC_by_morph_model.te <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + Iconicity.z*Grade + PSC_te.z*Grade*morph + (1|word), 
  data = df.phono
)
PSC_by_morph_model.te.aic = AIC(PSC_by_morph_model.te)
# 36742.81
summary(PSC_by_morph_model.te)
# target interaction is significant, but AIC isn't better
base.phono.aic - PSC_by_morph_model.te.aic
# -5.739636

quantile(df.phono$PSC_te.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
psc_by_morph_model.te.df = as.data.frame(
  ggpredict(
    PSC_by_morph_model.te, terms = c(
      "PSC_te.z[0.9284691, -0.5018487, -0.1012761, 0.3121661, 0.5781347, 0.7763081, 0.9571830, 1.1366814, 1.3214705]",
      "morph",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(psc_by_morph_model.te.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('PSC'[te])) +
  ylab('log(Freq)') +
  ggtitle(expression("3-way interaction of PSC[te], morphological complexity and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 



# Levenshtein neighbors
PSC_model.ld <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + PSC_ld.z*Grade + (1|word), 
  data = df.phono
)
PSC_model.ld.aic = AIC(PSC_model.ld)
# 35925.15
summary(PSC_model.ld)
# target interaction is significant
base.phono.aic - PSC_model.ld.aic
# 811.9283

quantile(df.phono$PSC_ld.z, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
PSC_model.ld.df = as.data.frame(
  ggpredict(
    PSC_model.ld, terms = c(
      "PSC_ld.z[-0.48016851, -0.05570006, 0.19447226, 0.37558923, 0.56505103, 0.80452613, 1.18087622, 1.65990238, 2.04141394]",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(PSC_model.ld.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('PSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("2-way interaction of PSC[ld] and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# PSCld*morph interaction
PSC_by_morph_model.ld <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + Iconicity.z*Grade + PSC_ld.z*Grade*morph + (1|word), 
  data = df.phono
)
PSC_by_morph_model.ld.aic = AIC(PSC_by_morph_model.ld)
# 35909.89
summary(PSC_by_morph_model.ld)
# target interaction is significant
PSC_model.ld.aic - PSC_by_morph_model.ld.aic
# 15.26154

PSC_by_morph_model.ld.df = as.data.frame(
  ggpredict(
    PSC_by_morph_model.ld, terms = c(
      "PSC_ld.z[-0.48016851, -0.05570006, 0.19447226, 0.37558923, 0.56505103, 0.80452613, 1.18087622, 1.65990238, 2.04141394]",
      "morph",
      "Grade[1,5,9,13]"
    )
  )
)

ggplot(PSC_by_morph_model.ld.df, aes(x = x, y = predicted, colour = as.factor(group), linetype = as.factor(group))) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(group)), alpha = 0.25) +
  xlab(expression('PSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("3-way interaction of PSC[ld], morphological complexity and Grade on Frequency")) +
  labs(colour = "Grade", fill = "Grade", linetype = 'Grade') +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

##### Interaction FSC*SND #####
# 2-way interaction: SND*OSC
osc_by_snd_model.ld <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_ld.z*Grade*snd.z + (1|word), 
  data = df.ortho
)

osc_by_snd_model.ld.aic = AIC(osc_by_snd_model.ld)
# 38341.46
summary(osc_by_snd_model.ld)
# target interaction is significant, but AIC only marginally better
osc_model.ld.aic - osc_by_snd_model.ld.aic
# 1.684696

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.ortho = data.frame(
  ggpredict(osc_by_snd_model.ld, terms = c("OSC_ld.z", "snd.z[-2.102646, -0.4752133, 1.1277678]", "Grade[1,5,9,13]"))
)
effects.2way_int.ortho$group <- mapvalues(
  effects.2way_int.ortho$group, 
  from = c("-2.102646", "-0.4752133", "1.1277678"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.ortho, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("3-way interaction of SND, OSC[ld] and Grade on log(Freq)")) +
  labs(colour = "SND", fill = "SND") +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# 2-way interaction: SND*PSC
psc_by_snd_model.ld <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + PSC_ld.z*Grade*snd.z + (1|word), 
  data = df.phono
)

psc_by_snd_model.ld.aic = AIC(psc_by_snd_model.ld)
# 35924.33
summary(psc_by_snd_model.ld)
# target interaction is significant, but AIC only slightly better
PSC_model.ld.aic - psc_by_snd_model.ld.aic
# 0.8138322

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.phono = data.frame(
  ggpredict(psc_by_snd_model.ld, terms = c("PSC_ld.z", "snd.z[-2.102646, -0.489336, 1.088724]", "Grade[1,5,9,13]"))
)
effects.2way_int.phono$group <- mapvalues(
  effects.2way_int.phono$group, 
  from = c("-2.102646", "-0.489336", "1.088724"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.phono, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('PSC'[ld])) +
  ylab('log(Freq)') +
  ggtitle(expression("3-way interaction of SND, PSC[ld] and Grade on log(Freq)")) +
  labs(colour = "SND", fill = "SND") +
  facet_grid(. ~ facet) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


##### RANDOM BASELINE #####
setwd("/media/gioca90/University/TiU/Research/Projects/ba_thesis/output/FSCrandom/")
iters = seq(1, 1000, 1)
stats = NULL
for (i in iters) {
  print(sprintf("Processing iteration %s", i))
  df_rnd = read.csv(sprintf("df%s.csv", i), header = T, sep = ',')  
  
  colnames(df_rnd)[which(names(df_rnd) == "OSC_te")] <- "OSC_te_rnd"
  colnames(df_rnd)[which(names(df_rnd) == "OSC_ld")] <- "OSC_ld_rnd"
  colnames(df_rnd)[which(names(df_rnd) == "PSC_te")] <- "PSC_te_rnd"
  colnames(df_rnd)[which(names(df_rnd) == "PSC_ld")] <- "PSC_ld_rnd"
  
  df.phono = merge(df.phono, df_rnd[c('word', "PSC_te_rnd", "PSC_ld_rnd")], by='word')
  df.ortho = merge(df.ortho, df_rnd[c('word', "OSC_te_rnd", "OSC_ld_rnd")], by='word')
  
  df.ortho$OSC_te_rnd.z = c(scale(box_cox_transf(df.ortho$OSC_te_rnd + 0.0001, -10, 10), center = T, scale = T))
  df.ortho$OSC_ld_rnd.z = c(scale(box_cox_transf(df.ortho$OSC_ld_rnd + 0.0001, -10, 10), center = T, scale = T))

  df.phono$PSC_te_rnd.z = c(scale(box_cox_transf(df.phono$PSC_te_rnd + 0.0001, -10, 10), center = T, scale = T))
  df.phono$PSC_ld_rnd.z = c(scale(box_cox_transf(df.phono$PSC_ld_rnd + 0.0001, -10, 10), center = T, scale = T))

  osc_model.te.rnd <- lmer(
    ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_te_rnd.z*Grade + (1|word), 
    data = df.ortho
    )
  deltaAIC.osc_te_rnd = base_model.ortho.int.aic - AIC(osc_model.te.rnd)
  
  osc_model.ld.rnd <- lmer(
    ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + old.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + OSC_ld_rnd.z*Grade + (1|word), 
    data = df.ortho
    )
  deltaAIC.osc_ld_rnd = base_model.ortho.int.aic - AIC(osc_model.ld.rnd)
  
  psc_model.te.rnd <- lmer(
    ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + PSC_te_rnd.z*Grade + (1|word), 
    data = df.phono
    )
  deltaAIC.psc_te_rnd = base_model.phono.int.aic - AIC(psc_model.te.rnd)
  
  psc_model.ld.rnd <- lmer(
    ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + pnd.z*snd.z*Grade + morph*Grade + Iconicity.z*Grade + PSC_ld_rnd.z*Grade + (1|word), 
    data = df.phono
    )
  deltaAIC.psc_ld_rnd = base_model.phono.int.aic - AIC(psc_model.ld.rnd)
  
  df.phono = subset(df.phono, select = -c(PSC_te_rnd, PSC_ld_rnd))
  df.ortho = subset(df.ortho, select = -c(OSC_te_rnd, OSC_ld_rnd))
  
  stats = rbind(stats,
                data.frame(i, deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd))
}

rm(i, deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd)

subtitle_osc_rnd.te = expression(OSC[te])
subtitle_osc_rnd.ld = expression(OSC[ld])
subtitle_psc_rnd.te = expression(PSC[te])
subtitle_psc_rnd.ld = expression(PSC[ld])

# DeltaAIC, linear models
range(stats$deltaAIC.osc_te_rnd)
p1.fsc_rnd.aic = ggplot(data = stats, aes(x=deltaAIC.osc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("log(Freq) and FSC: LMER " ~ Delta[AIC]), subtitle = subtitle_osc_rnd.te)  +
  labs(x = ' ', y = 'freq') +
  coord_cartesian(xlim = (c(-22, 77)))
(1+sum(stats$deltaAIC.osc_te_rnd >= (base_model.ortho.int.aic - osc_model.te.aic)))/(length(iters)+1)
# 0.006993007

range(stats$deltaAIC.osc_ld_rnd)
p2.fsc_rnd.aic = ggplot(data = stats, aes(x=deltaAIC.osc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_rnd.ld)  +
  labs(x = ' ', y = ' ')  +
  coord_cartesian(xlim = (c(-21, 1030)))
(1+sum(stats$deltaAIC.osc_ld_rnd >= (base_model.ortho.int.aic - osc_model.ld.aic)))/(length(iters)+1)
# 0.000999001

range(stats$deltaAIC.psc_te_rnd)
p3.fsc_rnd.aic = ggplot(data = stats, aes(x=deltaAIC.psc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - PSC_model.te.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq')  +
  coord_cartesian(xlim = (c(-22, 56)))
(1+sum(stats$deltaAIC.psc_te_rnd >= (base_model.phono.int.aic - PSC_model.te.aic)))/(length(iters)+1)
# 0.4105894

range(stats$deltaAIC.psc_ld_rnd)
p4.fsc_rnd.aic = ggplot(data = stats, aes(x=deltaAIC.psc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - PSC_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')   +
  coord_cartesian(xlim = (c(-21, 815)))
(1+sum(stats$deltaAIC.psc_ld_rnd >= (base_model.phono.int.aic - PSC_model.ld.aic)))/(length(iters)+1)
# 0.000999001

grid.arrange(p1.fsc_rnd.aic, p2.fsc_rnd.aic, p3.fsc_rnd.aic, p4.fsc_rnd.aic, ncol = 2)

range(stats$deltaAIC.psc_ld_rnd)

##### PCA FOR NEIGHBORHOOD MEASURES - ORTHOGRAPHY #####
neighborhoods_te.ortho.pca <- prcomp(df.ortho[c("OSC_te.z", "old.z", "snd.z")], center = TRUE, scale. = TRUE)
neighborhoods_ld.ortho.pca <- prcomp(df.ortho[c("OSC_ld.z", "old.z", "snd.z")], center = TRUE, scale. = TRUE)

plot(cumsum(neighborhoods_te.ortho.pca$sdev^2 / sum(neighborhoods_te.ortho.pca$sdev^2)), type="b", col='red')
lines(cumsum(neighborhoods_ld.ortho.pca$sdev^2 / sum(neighborhoods_ld.ortho.pca$sdev^2)), type="b", col='steelblue')

corrplot(cor(neighborhoods_te.ortho.pca$x[,1:3], df.ortho[c("OSC_te.z", "old.z", "snd.z")]))
corrplot(cor(neighborhoods_ld.ortho.pca$x[,1:3], df.ortho[c("OSC_ld.z", "old.z", "snd.z")]))

df.ortho['PC1_te_ortho'] = neighborhoods_te.ortho.pca$x[,1]
df.ortho['PC2_te_ortho'] = neighborhoods_te.ortho.pca$x[,2]
df.ortho['PC3_te_ortho'] = neighborhoods_te.ortho.pca$x[,3]

df.ortho['PC1_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,1]
df.ortho['PC2_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,2]
df.ortho['PC3_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,3]


pca_te_model.ortho <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + morph*Grade + PC1_te_ortho*Grade + PC2_te_ortho*Grade + PC3_te_ortho*Grade + (1|word), 
  data = df.ortho
)
pca_te_model.ortho.aic = AIC(pca_te_model.ortho)
# 39756.44
summary(pca_te_model.ortho)
# interactions are significant with all three principal components:
# PC1 high when OSC low, old low, and snd high; negative effect on Zeno freq at earliest grade, negative interaction with grade
# PC2 high when old high and snd high; negative effect on Zeno freq at earliest grade, negative interaction with grade
# PC3 high OSC low, old higher, snd lower; positive effect on Zeno freq at earliest grade, negative interaction with grade
base.ortho.aic - pca_te_model.ortho.aic
# -387.5999: AIC is however much worse than the base model


pca_ld_model.ortho <- lmer(
  ZenoFrequency ~ n_chars.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + morph*Grade + PC1_ld_ortho*Grade + PC2_ld_ortho*Grade + PC3_ld_ortho*Grade + (1|word), 
  data = df.ortho
)
pca_ld_model.ortho.aic = AIC(pca_ld_model.ortho)
# 38540.84
summary(pca_ld_model.ortho)
# interactions are significant with all three principal components:
# PC1 high when OSC high, OLD low, snd lower; positive effect on Zeno freq at earliest grade, negative interaction with grade
# PC2 high when OLD higher, snd low; positive effect on Zeno freq at earliest grade, negative interaction with grade
# PC3 high when OSC lower, old lower, snd slightly lower; no effect on Zeno freq, positive interaction
base.ortho.aic - pca_ld_model.ortho.aic
# 828.0069


##### PCA FOR NEIGHBORHOOD MEASURES - PHONOLOGY #####
neighborhoods_te.phono.pca <- prcomp(df.phono[c("PSC_te.z", "pnd.z", "snd.z")], center = TRUE, scale. = TRUE)
neighborhoods_ld.phono.pca <- prcomp(df.phono[c("PSC_ld.z", "pnd.z", "snd.z")], center = TRUE, scale. = TRUE)

plot(cumsum(neighborhoods_te.phono.pca$sdev^2 / sum(neighborhoods_te.phono.pca$sdev^2)), type="b", col='red')
lines(cumsum(neighborhoods_ld.phono.pca$sdev^2 / sum(neighborhoods_ld.phono.pca$sdev^2)), type="b", col='steelblue')

corrplot(cor(neighborhoods_te.phono.pca$x[,1:3], df.phono[c("PSC_te.z", "pnd.z", "snd.z")]))
corrplot(cor(neighborhoods_ld.phono.pca$x[,1:3], df.phono[c("PSC_ld.z", "pnd.z", "snd.z")]))

df.phono['PC1_te_phono'] = neighborhoods_te.phono.pca$x[,1]
df.phono['PC2_te_phono'] = neighborhoods_te.phono.pca$x[,2]
df.phono['PC3_te_phono'] = neighborhoods_te.phono.pca$x[,3]

df.phono['PC1_ld_phono'] = neighborhoods_ld.phono.pca$x[,1]
df.phono['PC2_ld_phono'] = neighborhoods_ld.phono.pca$x[,2]
df.phono['PC3_ld_phono'] = neighborhoods_ld.phono.pca$x[,3]


pca_te_model.phono <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + morph*Grade + PC1_te_phono*Grade + PC2_te_phono*Grade + PC3_te_phono*Grade + (1|word), 
  data = df.phono
)
pca_te_model.phono.aic = AIC(pca_te_model.phono)
# 37217.85
summary(pca_te_model.phono)
# interactions with grade are significant for all three principal components:
# PC1 high when PSC low, pnd lower, snd high: negative effect on Zeno freq at earliest grade, positive interaction with grade
# PC2 high when PSC high, pnd low: negative effect on Zeno freq at earliest grade, positive interaction with grade
# PC3 high when PSC low, pnd low, snd low: positive effect on Zeno freq, negative interaction with grade
base.phono.aic - pca_te_model.phono.aic
# -480.7717: AIC much lower than base model


pca_ld_model.phono <- lmer(
  ZenoFrequency ~ n_phons.z*Grade + concr.z*Grade + val.z*Grade + freq.z*Grade + morph*Grade + PC1_ld_phono*Grade + PC2_ld_phono*Grade + PC3_ld_phono*Grade + (1|word), 
  data = df.phono
)
pca_ld_model.phono.aic = AIC(pca_ld_model.phono)
# 36229.56
summary(pca_ld_model.phono)
# interactions with grade are significant for all three principal components:
# PC1 high when PSC high, pnd high, snd low: positive effect on Zeno freq at earliest grade, negative interaction with grade
# PC2 high when pnd lower and snd low: positive effect on Zeno freq at earliest grade, negative interaction with grade
# PC3 high when PSC higher, pnd lower: negative effect on Zeno freq at earliest grade, positive interaction with grade
base.phono.aic - pca_ld_model.phono.aic
# 507.514
