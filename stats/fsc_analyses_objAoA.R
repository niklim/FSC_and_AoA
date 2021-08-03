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
library(readxl)
library(mgcv)
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
fsc_data <- read.csv(file = "output/fsc_measures.csv")
iconicity <- read.csv(file = '/media/gioca90/University/TiU/Research/Resources/iconicity_ratings.csv')
df.mald = read.csv(file = '/media/gioca90/University/TiU/Research/Resources/MALD/AllData.txt', sep = '\t')
bb_aoa = read_excel('data/Master file with all values for test based AoA measures.xlsx')

data = merge(fsc_data, bb_aoa[c('WORD', 'AoAtestbased')], by.x='word', by.y='WORD')
data_sorted <- data[order(data$word, data$AoAtestbased, decreasing = FALSE),]
data_unique <- data_sorted[ !duplicated(data_sorted$word), ] 

#transform morph into factor
data_unique$morph <- factor(data_unique$morph,   
                            level=c(0,1),    
                            labels=c("Mono","Poly"))


### Transform variables into normal shape and comparable units (SD) ###
df.ortho <- data_unique[,FALSE]
df.phono <- data_unique[,FALSE]

df.ortho$word = data_unique$word
df.phono$word = data_unique$word

df.ortho$AoA = data_unique$AoAtestbased
df.phono$AoA = data_unique$AoAtestbased

df.ortho$morph = data_unique$morph
df.phono$morph = data_unique$morph

df.ortho$concr.z = c(scale(box_cox_transf(data_unique$concr, -10, 10), center = T, scale = T))
df.phono$concr.z = df.ortho$concr.z

df.ortho$val.z = c(scale(box_cox_transf(data_unique$val, -10, 10), center = T, scale = T))
df.phono$val.z = df.ortho$val.z

df.ortho$freq.z = c(scale(box_cox_transf(data_unique$freq, -10, 10), center = T, scale = T))
df.phono$freq.z = df.ortho$freq.z

df.ortho$snd.z = c(scale(box_cox_transf(data_unique$snd, -10, 10), center = T, scale = T))
df.phono$snd.z = df.ortho$snd.z

df.ortho$n_chars.z = c(scale(box_cox_transf(data_unique$length, -10, 10), center = T, scale = T))
df.ortho$old.z = c(scale(box_cox_transf(data_unique$old20, -10, 10), center = T, scale = T))
df.ortho$OSC_te.z = c(scale(box_cox_transf(data_unique$OSC_te + 0.0001, -10, 10), center = T, scale = T))
df.ortho$OSC_ld.z = c(scale(box_cox_transf(data_unique$OSC_ld, -10, 10), center = T, scale = T))

df.phono$n_phons.z = c(scale(box_cox_transf(data_unique$n_phon, -10, 10), center = T, scale = T))
df.phono$PSC_te.z = c(scale(box_cox_transf(data_unique$PSC_te + 0.0001, -10, 10), center = T, scale = T))
df.phono$PSC_ld.z = c(scale(box_cox_transf(data_unique$PSC_ld, -10, 10), center = T, scale = T))
df.phono = merge(df.phono, df.mald[c('Item', 'PhonND')], by.x='word', by.y='Item')
df.phono = distinct(df.phono)
df.phono$pnd.z = c(scale(box_cox_transf(df.phono$PhonND + 0.0001, -10, 10), center = T, scale = T))

df.ortho.icon = merge(df.ortho, iconicity, by.x='word', by.y='Word')
df.ortho.icon$iconicity.z = c(
  scale(box_cox_transf(df.ortho.icon$Iconicity + abs(min(df.ortho.icon$Iconicity)) + 0.0001, -10, 10), 
        center = T, scale = T)
)

df.phono.icon = merge(df.phono, iconicity, by.x='word', by.y='Word')
df.phono.icon$iconicity.z = c(
  scale(box_cox_transf(df.phono.icon$Iconicity + abs(min(df.phono.icon$Iconicity)) + 0.0001, -10, 10), 
        center = T, scale = T)
)

rm(bb_aoa, data, data_sorted, data_unique, df.mald, fsc_data, iconicity)

phono <- df.phono[c("AoA", "PSC_te.z", "PSC_ld.z", "concr.z", "val.z", "n_phons.z", "freq.z", "pnd.z", "snd.z")]
colnames(phono) <- c("AoA", "PSCte", "PSCld", "Concr", "Val", "Len", "Freq", "PND", "SND")
pairs.panels(phono,
             method = "pearson",
             hist.col = "steelblue",
             density = TRUE,
             ellipses = FALSE,
             pch = '.',
             cex = 2.5,
             cex.labels = 2,
             cex.main = 2,
             main = 'Phonological variables - Objective AoA')
rm(phono)

ortho <- df.ortho[c("AoA", "OSC_te.z", "OSC_ld.z", "concr.z", "val.z", "n_chars.z", "freq.z", "old.z", "snd.z")]
colnames(ortho) <- c("AoA", "OSCte", "OSCld", "Concr", "Val", "Len", "Freq", "OLD20", "SND")
pairs.panels(ortho,
             method = "pearson",
             hist.col = "steelblue",
             density = TRUE,
             ellipses = FALSE,
             pch = '.',
             cex = 2.5,
             cex.labels = 2,
             cex.main = 2,
             main = "Orthographic variables - Objective AoA")
rm(ortho)

##### MULTIPLE REGRESSION - ORTHOGRAPHY #####
# Base Model
base_model.ortho <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph, data = df.ortho
)
summary(base_model.ortho)$r.sq
# 0.4347285
summary(base_model.ortho)
# all predictors are significant
VIF(base_model.ortho)
# some collinearity between length and OLD
base.ortho.aic <- AIC(base_model.ortho)
# 38433.19

# Base Model with snd*old interaction
base_model.ortho.int <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z * snd.z, data = df.ortho
)
base_model.ortho.int.aic = AIC(base_model.ortho.int)
# 38434.26
summary(base_model.ortho.int)$r.sq
# 0.4347948
summary(base_model.ortho.int)
# the target interaction is not significant
base.ortho.aic - base_model.ortho.int.aic
# -1.068727

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
quantile(df.ortho$old.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
base_model.ortho.int.df = as.data.frame(
  ggpredict(
    base_model.ortho.int, terms = c(
      "old.z[-1.8499734, -1.4825293, -0.7839768, -0.4549141, -0.3046405, -0.1624940, 0.2820197, 0.6621340, 0.8569692, 1.3284692, 3.5165398]",
      "snd.z[-1.6252009, -0.02585081, 1.65709681]")
  )
)
base_model.ortho.int.df$group <- mapvalues(
  base_model.ortho.int.df$group, 
  from = c("-1.6252009", "-0.02585081", "1.65709681"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.ortho.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('OLD') +
  ylab('AoA') +
  ggtitle(expression("Interaction of OLD20 and SND on objective AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Iconicity control
icon_model.ortho <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + iconicity.z, data = df.ortho.icon
)
icon_model.ortho.aic = AIC(icon_model.ortho)
# 8336.735
summary(icon_model.ortho)$r.sq
# 0.4793029
summary(icon_model.ortho)
# iconicity is significant, snd and length aren't
VIF(icon_model.ortho)
# no worrying effects of collinearity


# target embedded neighbors
osc_model.te <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + OSC_te.z, data = df.ortho
)
osc_model.te.aic = AIC(osc_model.te)
# 38415.14
summary(osc_model.te)$r.sq
# 0.4361561
summary(osc_model.te)
# OSCte is significant
base.ortho.aic - osc_model.te.aic
# 18.05256


# target embedded neighbors over iconicity
osc_icon_model.te <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + iconicity.z + OSC_te.z, data = df.ortho.icon
)
osc_icon_model.te.aic = AIC(osc_icon_model.te)
# 8307.486
summary(osc_icon_model.te)$r.sq
# 0.4873794
summary(osc_icon_model.te)
# both iconicity and OSCte are significant
icon_model.ortho.aic - osc_icon_model.te.aic
# 29.24947


# OSCte*morph interaction
osc_by_morph_model.te <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph*OSC_te.z, data = df.ortho
)
osc_by_morph_model.te.aic = AIC(osc_by_morph_model.te)
# 38366.12
summary(osc_by_morph_model.te)$r.sq
# 0.4397717
summary(osc_by_morph_model.te)
# target interaction is significant
osc_model.te.aic - osc_by_morph_model.te.aic
# 49.01498

quantile(df.ortho$OSC_te.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
osc_by_morph_model.te.df = as.data.frame(
  ggpredict(
    osc_by_morph_model.te, terms = c(
      "OSC_te.z[-2.2842981, -1.4623802, -1.1448036, -0.6867680, -0.1620119, 0.2259678, 0.4907714, 0.7114467, 0.9463726, 1.2089950, 1.9716403]",
      "morph")
  )
)

ggplot(osc_by_morph_model.te.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('OSC'[te])) +
  ylab('AoA') +
  ggtitle(expression("Interaction of" ~ OSC[te] ~ "and morphological status on AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Levenshtein neighbors
osc_model.ld <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + OSC_ld.z, data = df.ortho
)
osc_model.ld.aic = AIC(osc_model.ld)
# 38302.7
summary(osc_model.ld)$r.sq
# 0.4440944
summary(osc_model.ld)
# OSCld is a significant predictor
base.ortho.aic - osc_model.ld.aic
# 130.4921


# Levenshtein neighbors over iconicity
osc_icon_model.ld <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + iconicity.z + OSC_ld.z, data = df.ortho.icon
)
osc_icon_model.ld.aic = AIC(osc_icon_model.ld)
# 8320.334
summary(osc_icon_model.ld)$r.sq
# 0.4840741
summary(osc_icon_model.ld)
# both iconicity and OSCld are significant predictors
icon_model.ortho.aic - osc_icon_model.ld.aic
# 16.40152


# OSCld*morph interaction
osc_by_morph_model.ld <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph*OSC_ld.z, data = df.ortho
)
osc_by_morph_model.ld.aic = AIC(osc_by_morph_model.ld)
# 38296.87
summary(osc_by_morph_model.ld)$r.sq
# 0.444643
summary(osc_by_morph_model.ld)
# interaction is significant
osc_model.ld.aic - osc_by_morph_model.ld.aic
# 5.82956

quantile(df.ortho$OSC_ld.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
osc_by_morph_model.ld.df = as.data.frame(
  ggpredict(
    osc_by_morph_model.ld, terms = c(
      "OSC_ld.z[-3.5470053, -1.3026135, -0.8370214, -0.5002286, -0.2273463, 0.0241192, 0.2703450, 0.5157738, 0.8055791, 1.2555764, 3.1631114]",
      "morph")
  )
)

ggplot(osc_by_morph_model.ld.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('AoA') +
  ggtitle(expression("Interaction of" ~ OSC[ld] ~ "and morphological status on AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


##### MULTIPLE REGRESSION - PHONOLOGY #####
# Base Model
base_model.phono <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph, data = df.phono
)
summary(base_model.phono)$r.sq
# 0.4580349
summary(base_model.phono)
# all predictors are significant
VIF(base_model.phono)
# only moderate collinearity between number of phonemes and PND
base.phono.aic <- AIC(base_model.phono)
# 30442.48


# Base Model with pnd*snd interaction
base_model.phono.int <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph, data = df.phono
)
summary(base_model.phono.int)$r.sq
# 0.45815
summary(base_model.phono.int)
# target interaction isn't significant
VIF(base_model.phono.int)
# no extra collinearity risk
base_model.phono.int.aic <- AIC(base_model.phono.int)
# 30443.12
base.phono.aic - base_model.phono.int.aic
# -0.6394581

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
quantile(df.phono$pnd.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
base_model.phono.int.df = as.data.frame(
  ggpredict(
    base_model.phono.int, terms = c(
      "pnd.z[-2.16071576, -0.85287765, -0.85287765, -0.47100090, -0.35634542, -0.05141835, 0.29885338, 0.61930081, 0.93519895, 1.30749626, 2.16767218 ]",
      "snd.z[-1.65190024, -0.04120388, 1.60161652]")
  )
)
base_model.phono.int.df$group <- mapvalues(
  base_model.phono.int.df$group, 
  from = c("-1.65190024", "-0.04120388", "1.60161652"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.phono.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('PND') +
  ylab('AoA') +
  ggtitle(expression("Interaction of PND and SND on objective AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# Iconicity control
icon_model.phono <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + iconicity.z, data = df.phono.icon
)
icon_model.phono.aic = AIC(icon_model.phono)
# 7229.584
summary(icon_model.phono)$r.sq
# 0.4760113
summary(icon_model.phono)
# iconicity is a significant predictor
cor.test(df.phono.icon$iconicity.z, df.phono.icon$AoA)

# target embedded neighbors
psc_model.te <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + PSC_te.z, data = df.phono
)
psc_model.te.aic = AIC(psc_model.te)
# 30426.31
summary(psc_model.te)$r.sq
# 0.4595699
summary(psc_model.te)
# PSCte is a significant predictor
base.phono.aic - psc_model.te.aic
# 16.17232


# target embedded neighbors over iconicity
psc_icon_model.te <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + iconicity.z + PSC_te.z, data = df.phono.icon
)
psc_icon_model.te.aic = AIC(psc_icon_model.te)
# 7227.144
summary(psc_icon_model.te)$r.sq
# 0.4773233
summary(psc_icon_model.te)
# iconicity is significant, PSCte only slightly
icon_model.phono.aic - psc_icon_model.te.aic
# 2.43962

# PSCte*morph interaction
psc_by_morph_model.te <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph*PSC_te.z, data = df.phono
)
psc_by_morph_model.te.aic = AIC(psc_by_morph_model.te)
# 30404.1
summary(psc_by_morph_model.te)$r.sq
# 0.4616077
summary(psc_by_morph_model.te)
# target interaction is significant
psc_model.te.aic - psc_by_morph_model.te.aic
# 22.20384

quantile(df.phono$PSC_te.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
psc_by_morph_model.te.df = as.data.frame(
  ggpredict(
    psc_by_morph_model.te, terms = c(
      "PSC_te.z[-2.1897122, -1.3872573, -1.1133175, -0.7636885, -0.2131268, 0.2770647, 0.5643153, 0.7923830, 1.0100233, 1.2435375, 1.8546669]",
      "morph")
  )
)

ggplot(psc_by_morph_model.te.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('PSC'[te])) +
  ylab('AoA') +
  ggtitle(expression("Interaction of" ~ PSC[te] ~ "and morphological status on AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Levenshtein neighbors
psc_model.ld <- lm(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + PSC_ld.z, data = df.phono)
psc_model.ld.aic = AIC(psc_model.ld)
# 30194.1
summary(psc_model.ld)$r.sq
# 0.4788061
summary(psc_model.ld)
# PSCld is significant
base.phono.aic - psc_model.ld.aic
# 248.3814

# Levenshtein neighbors over iconicity
psc_icon_model.ld <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + iconicity.z + PSC_ld.z, data = df.phono.icon
)
psc_icon_model.ld.aic = AIC(psc_icon_model.ld)
# 7213.149
summary(psc_icon_model.ld)$r.sq
# 0.4814375
summary(psc_icon_model.ld)
# both iconicity and PSCld are significant
icon_model.phono.aic - psc_icon_model.ld.aic
# 16.4352

# PSCld*morph interaction
psc_by_morph_model.ld <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph*PSC_ld.z, data = df.phono
)
psc_by_morph_model.ld.aic = AIC(psc_by_morph_model.ld)
# 30195.96
summary(psc_by_morph_model.ld)$r.sq
# 0.4788168
summary(psc_by_morph_model.ld)
# target interaction is not significant
psc_model.ld.aic - psc_by_morph_model.ld.aic
# -1.867567

##### Interaction FSC*SND #####
### replication iconicity*SND interaction with FSC measures ###

cor.test(df.ortho.icon$iconicity.z, df.ortho.icon$snd.z)
# 0.06485585 [0.02107388 0.10838956], t = 2.9044, df = 1997, p-value = 0.00372
cor.test(df.phono.icon$iconicity.z, df.phono.icon$snd.z)
# 0.06377966 [0.01725166 0.11003204], t = 2.688, df = 1769, p-value = 0.007255

cor.test(df.ortho.icon$iconicity.z, df.ortho.icon$OSC_ld.z)
# 0.2074537 [0.1651136 0.2490307], t = 9.4768, df = 1997, p-value < 2.2e-16
cor.test(df.phono.icon$iconicity.z, df.phono.icon$PSC_ld.z)
# 0.1610802 [0.1153665 0.2061129], t = 6.8646, df = 1769, p-value = 9.196e-12

cor.test(df.ortho$OSC_ld.z, df.ortho$snd.z)
# -0.2885349 [-0.3085853 -0.2682281], t = -26.832, df = 7928, p-value < 2.2e-16
cor.test(df.phono$PSC_ld.z, df.phono$snd.z)
# -0.3259303 [-0.3476427 -0.3038684], t = -27.591, df = 6405, p-value < 2.2e-16

# predict systematicity (FSC-ld)
OSC_ld.lm <- lm(formula = OSC_ld.z ~ n_chars.z + freq.z + old.z + snd.z + AoA, data = df.ortho)
summary(OSC_ld.lm)

PSC_ld.lm <- lm(formula = PSC_ld.z ~ n_phons.z + freq.z + pnd.z + snd.z + AoA, data = df.phono)
summary(PSC_ld.lm)


# Interaction: SND*OSC
osc_by_snd_model.ld <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + OSC_ld.z * snd.z, data = df.ortho
)
osc_by_snd_model.ld.aic = AIC(osc_by_snd_model.ld)
# 38291.14
summary(osc_by_snd_model.ld)$r.sq
# 0.4450438
summary(osc_by_snd_model.ld)
# target interaction is significant
osc_model.ld.aic - osc_by_snd_model.ld.aic
# 11.55418

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.ortho = data.frame(
  ggpredict(osc_by_snd_model.ld, terms = c("OSC_ld.z","snd.z[-1.6252009, -0.02585081, 1.65709681]"))
)
effects.2way_int.ortho$group <- mapvalues(
  effects.2way_int.ortho$group, 
  from = c("-1.6252009", "-0.02585081", "1.65709681"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.ortho, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('AoA') +
  ggtitle(expression("Interaction of SND and" ~ OSC[ld] ~ "on AoA")) +
  labs(colour = "SND", fill = "SND")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# Interaction: SND*PSC
psc_by_snd_model.ld <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + PSC_ld.z * snd.z, data = df.phono
)
psc_by_snd_model.ld.aic = AIC(psc_by_snd_model.ld)
# 30163.76
summary(psc_by_snd_model.ld)$r.sq
# 0.4814298
summary(psc_by_snd_model.ld)
# target interaction is significant
psc_model.ld.aic - psc_by_snd_model.ld.aic
# 30.33438

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.phono = data.frame(
  ggpredict(psc_by_snd_model.ld, terms = c("PSC_ld.z","snd.z[-1.65190024, -0.04120388, 1.60161652]"))
)
effects.2way_int.phono$group <- mapvalues(
  effects.2way_int.phono$group, 
  from = c("-1.65190024", "-0.04120388", "1.60161652"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.phono, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('PSC'[ld])) +
  ylab('AoA') +
  ggtitle(expression("Interaction of SND and" ~ PSC[ld] ~ "on AoA")) +
  labs(colour = "SND", fill = "SND")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


##### GAMs #####
# Base gam
base_gam.phono <- bam(AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph, data = df.phono)
summary(base_gam.phono)
# all smooths are significant
base.phono.aic.gam <- AIC(base_gam.phono)

par(mfrow=c(2,3))
plot_smooth(base_gam.phono, view="freq.z", rug=TRUE, col="steelblue", main = "Phonology predictors (base gam)",
            xlab = 'frequency', ylab = 'AoA')
plot_smooth(base_gam.phono, view="pnd.z", rug=TRUE, col="steelblue", 
            xlab = 'PND', ylab = 'AoA')
plot_smooth(base_gam.phono, view="n_phons.z", rug=TRUE, col="steelblue", 
            xlab = 'Length', ylab = 'AoA')
plot_smooth(base_gam.phono, view="concr.z", rug=TRUE, col="steelblue", 
            xlab = 'Concreteness', ylab = 'AoA')
plot_smooth(base_gam.phono, view="val.z", rug=TRUE, col="steelblue", 
            xlab = 'Valence', ylab = 'AoA')
plot_smooth(base_gam.phono, view="snd.z", rug=TRUE, col="steelblue", 
            xlab = 'SND', ylab = 'AoA')


# Base gam with pnd*snd interaction
base_gam.phono.int <- bam(AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + ti(pnd.z, snd.z) + morph, 
                          data = df.phono)
summary(base_gam.phono.int)
# tensor product not significant
base_gam.phono.int.aic <- AIC(base_gam.phono.int)
base.phono.aic.gam - base_gam.phono.int.aic
# -0.5016009

par(mfrow=c(1, 1))
vis.gam(base_gam.phono.int, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'pnd.z'), main = "composite smooth: PND, SND", xlab = "SND", ylab = "PND")


# Iconicity control
icon_gam.phono <- bam(
  AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z), data = df.phono.icon
)
icon_gam.phono.aic = AIC(icon_gam.phono)
summary(icon_gam.phono)
# iconicity smooth is significant

plot_smooth(icon_gam.phono, view="iconicity.z", rug=TRUE, col="steelblue", main='Effect of iconicity on AoA (phono)',
            xlab = 'Iconicity', ylab = 'AoA')


# target embedded neighbors
psc_gam.te <- bam(AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(PSC_te.z), data = df.phono)
psc_gam.te.aic = AIC(psc_gam.te)
summary(psc_gam.te)
# smooth for PSC_te is significant
base.phono.aic.gam - psc_gam.te.aic
# 46.48438

plot_smooth(psc_gam.te, view="PSC_te.z", rug=TRUE, col="steelblue", main='Effect of PSC[te] on AoA', xlab = expression('PSC'[te]), ylab = 'AoA')


# target embedded neighbors over iconicity
psc_icon_gam.te <- bam(
  AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z) + s(PSC_te.z), 
  data = df.phono.icon
)
psc_icon_gam.te.aic = AIC(psc_icon_gam.te)
summary(psc_icon_gam.te)
# smooth for PSCte significant after including iconicity
icon_gam.phono.aic - psc_icon_gam.te.aic
# 6.103832

par(mfrow=c(1, 2))
plot_smooth(psc_icon_gam.te, view="PSC_te.z", rug=TRUE, col="steelblue", 
            main='Effect of PSC[te] abd Iconicity on AoA', xlab = expression('PSC'[te]), ylab = 'AoA')
plot_smooth(psc_icon_gam.te, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')



# Levenshtein neighbors
psc_gam.ld <- bam(AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(PSC_ld.z), data = df.phono)
psc_gam.ld.aic = AIC(psc_gam.ld)
summary(psc_gam.ld)
# smooth for PSCld is significant
base.phono.aic.gam - psc_gam.ld.aic
# 285.0209

par(mfrow=c(1, 1))
plot_smooth(psc_gam.ld, view="PSC_ld.z", rug=TRUE, col="steelblue", main='Effect of PSC[ld] on AoA', xlab = expression('PSC'[ld]), ylab = 'AoA')


# Levenshtein neighbors over iconicity
psc_icon_gam.ld <- bam(
  AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z) + s(PSC_ld.z), 
  data = df.phono.icon
)
psc_icon_gam.ld.aic = AIC(psc_icon_gam.ld)
summary(psc_icon_gam.ld)
# both smooths are significant
icon_gam.phono.aic - psc_icon_gam.ld.aic
# 38.45144

par(mfrow=c(1, 2))
plot_smooth(psc_icon_gam.ld, view="PSC_ld.z", rug=TRUE, col="steelblue", 
            main='Effect of PSC[ld] abd Iconicity on AoA', xlab = expression('PSC'[ld]), ylab = 'AoA')
plot_smooth(psc_icon_gam.ld, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')


# Interaction: SND*PSC
psc_by_snd_gam.ld <- bam(
    AoA ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + morph + snd.z*PSC_ld.z, 
    data = df.phono
)

summary(psc_by_snd_gam.ld)
# not significant
psc_by_snd_gam.ld.aic <- AIC(psc_by_snd_gam.ld)
psc_gam.ld.aic - psc_by_snd_gam.ld.aic
# -1.918932

par(mfrow=c(1, 1))
vis.gam(base_gam.phono.int, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'pnd.z'), main = "composite smooth: PND, SND", xlab = "SND", ylab = "PND")

##### Random Forest Regression Models #####

### Random Forest Base ###
rf_base.ortho <- ranger(formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z + snd.z, data = df.ortho, mtry = 5/3, importance = "impurity" )
print(rf_base.ortho)
vip(rf_base.ortho, aesthetics = list(col = "darkblue", fill = "darkblue"))
rf_base.ortho.rsq = rf_base.ortho$r.squared

rf_base.phono <- ranger(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + morph + pnd.z + snd.z, data = df.phono, mtry = 5/3, importance = "impurity" )
print(rf_base.phono)
vip(rf_base.phono, aesthetics = list(col = "darkblue", fill = "darkblue"))
rf_base.phono.rsq = rf_base.phono$r.squared

subtitle_osc.ld = expression(OSC[ld])
subtitle_osc.te = expression(OSC[te])
subtitle_psc.ld = expression(PSC[ld])
subtitle_psc.te = expression(PSC[te])


### Random Forest OSCte ###
rf_OSC_te <- ranger(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + OSC_te.z, data = df.ortho, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_te)
names(rf_OSC_te$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "OSC[te]"
))
rf_OSC_te.vip_df = data.frame(
  var=names(rf_OSC_te$variable.importance), Importance=rf_OSC_te$variable.importance, row.names=NULL
)
subtitle_osc.te = parse(text = "base[ortho] + OSC[te]")
p1.rf = ggplot(data = rf_OSC_te.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs) - objective AoA", subtitle = subtitle_osc.te) +
  coord_flip(ylim = c(0, 20000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_OSC_te.vip_df$var, rf_OSC_te.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_OSC_te$r.squared - rf_base.ortho.rsq


### Random Forest OSCld ###
rf_OSC_ld <- ranger(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + OSC_ld.z, data = df.ortho, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_ld)
names(rf_OSC_ld$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "OSC[ld]"
))
rf_OSC_ld.vip_df = data.frame(
  var=names(rf_OSC_ld$variable.importance), Importance=rf_OSC_ld$variable.importance, row.names=NULL
)
subtitle_osc.ld = parse(text = "base[ortho] + OSC[ld]")
p2.rf = ggplot(data = rf_OSC_ld.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle(" ", subtitle = subtitle_osc.ld) +
  coord_flip(ylim = c(0, 20000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_OSC_ld.vip_df$var, rf_OSC_ld.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_OSC_ld$r.squared - rf_base.ortho.rsq


### Random Forest PSCte ###
rf_PSC_te <- ranger(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + PSC_te.z, data = df.phono, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_te)
names(rf_PSC_te$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "PSC[te]"
))
rf_PSC_te.vip_df = data.frame(
  var=names(rf_PSC_te$variable.importance), Importance=rf_PSC_te$variable.importance, row.names=NULL
)
subtitle_psc.te = parse(text = "base[phono] + PSC[te]")
p3.rf = ggplot(data = rf_PSC_te.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs) - objective AoA", subtitle = subtitle_psc.te) +
  coord_flip(ylim = c(0, 20000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_PSC_te.vip_df$var, rf_PSC_te.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_PSC_te$r.squared - rf_base.phono.rsq


### Random Forest PSCld ###
rf_PSC_ld <- ranger(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + PSC_ld.z, data = df.phono, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_ld)
names(rf_PSC_ld$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "PSC[ld]"
))
rf_PSC_ld.vip_df = data.frame(
  var=names(rf_PSC_ld$variable.importance), Importance=rf_PSC_ld$variable.importance, row.names=NULL
)
subtitle_psc.ld = parse(text = "base[phono] + PSC[ld]")
p4.rf = ggplot(data = rf_PSC_ld.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("", subtitle = subtitle_osc.te) +
  coord_flip(ylim = c(0, 20000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_PSC_ld.vip_df$var, rf_PSC_ld.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_PSC_ld$r.squared - rf_base.phono.rsq

grid.arrange(p1.rf, p2.rf, ncol = 1)
grid.arrange(p3.rf, p4.rf, ncol = 1)



### Random Forest Iconicity Base ###
rf_base.ortho.icon <- ranger(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z + snd.z + iconicity.z, data = df.ortho.icon, mtry = 5/3, importance = "impurity"
)
print(rf_base.ortho.icon)
rf_base.ortho.icon.rsq = rf_base.ortho.icon$r.squared

rf_base.phono.icon <- ranger(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + morph + pnd.z + snd.z + iconicity.z, data = df.phono.icon, mtry = 5/3, importance = "impurity"
)
print(rf_base.phono.icon)
rf_base.phono.icon.rsq = rf_base.phono.icon$r.squared


subtitle_osc.ld.icon = expression("icon +" ~ OSC[ld])
subtitle_osc.te.icon = expression("icon +" ~ OSC[te])
subtitle_psc.ld.icon = expression("icon +" ~ PSC[ld])
subtitle_psc.te.icon = expression("icon +" ~ PSC[te])


### Random Forest Iconicity + OSCte ###
rf_OSC_te.icon <- ranger(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + iconicity.z + OSC_te.z, 
  data = df.ortho.icon, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_te.icon)
names(rf_OSC_te.icon$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "Icon", "OSC[te]"
))
rf_OSC_te.icon.vip_df = data.frame(
  var=names(rf_OSC_te.icon$variable.importance), Importance=rf_OSC_te.icon$variable.importance, row.names=NULL
)
subtitle_osc.te.icon = parse(text = "icon[ortho] + OSC[te]")
p1.rf.icon = ggplot(data = rf_OSC_te.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs w/ Iconicity) - objective AoA", subtitle = subtitle_osc.te.icon) +
  coord_flip(ylim = c(0, 3000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_OSC_te.icon.vip_df$var, rf_OSC_te.icon.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_OSC_te.icon$r.squared - rf_base.ortho.icon.rsq


### Random Forest Iconicity + OSCld ###
rf_OSC_ld.icon <- ranger(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + iconicity.z + OSC_ld.z, 
  data = df.ortho.icon, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_ld.icon)
names(rf_OSC_ld.icon$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "Icon", "OSC[ld]"
))
rf_OSC_ld.icon.vip_df = data.frame(
  var=names(rf_OSC_ld.icon$variable.importance), Importance=rf_OSC_ld.icon$variable.importance, row.names=NULL
)
subtitle_osc.ld.icon = parse(text = "icon[ortho] + OSC[ld]")
p2.rf.icon = ggplot(data = rf_OSC_ld.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle(" ", subtitle = subtitle_osc.ld.icon) +
  coord_flip(ylim = c(0, 3000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_OSC_ld.icon.vip_df$var, rf_OSC_ld.icon.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_OSC_ld.icon$r.squared - rf_base.ortho.icon.rsq

### Random Forest Iconicity + PSCte ###
rf_PSC_te.icon <- ranger(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + iconicity.z + PSC_te.z, 
  data = df.phono.icon, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_te.icon)
names(rf_PSC_te.icon$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "Icon", "PSC[te]"
))
rf_PSC_te.icon.vip_df = data.frame(
  var=names(rf_PSC_te.icon$variable.importance), Importance=rf_PSC_te.icon$variable.importance, row.names=NULL
)
subtitle_psc.te.icon = parse(text = "icon[phono] + PSC[te]")
p3.rf.icon = ggplot(data = rf_PSC_te.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs w/ Iconicity) - objective AoA", subtitle = subtitle_psc.te.icon) +
  coord_flip(ylim = c(0, 3000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_PSC_te.icon.vip_df$var, rf_PSC_te.icon.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_PSC_te.icon$r.squared - rf_base.phono.icon.rsq


### Random Forest Iconicity + PSCld ###
rf_PSC_ld.icon <- ranger(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + iconicity.z + PSC_ld.z, 
  data = df.phono.icon, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_ld.icon)
names(rf_PSC_ld.icon$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "Icon", "PSC[ld]"
))
rf_PSC_ld.icon.vip_df = data.frame(
  var=names(rf_PSC_ld.icon$variable.importance), Importance=rf_PSC_ld.icon$variable.importance, row.names=NULL
)
subtitle_psc.ld.icon = parse(text = "icon[phono] + PSC[ld]")
p4.rf.icon = ggplot(data = rf_PSC_ld.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle(" ", subtitle = subtitle_psc.ld.icon) +
  coord_flip(ylim = c(0, 3000)) +
  scale_x_discrete(
    labels = parse(text = levels(reorder(rf_PSC_ld.icon.vip_df$var, rf_PSC_ld.icon.vip_df$Importance)))
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
rf_PSC_ld.icon$r.squared - rf_base.phono.icon.rsq

grid.arrange(p1.rf.icon, p2.rf.icon, ncol = 1)
grid.arrange(p3.rf.icon, p4.rf.icon, ncol = 1)


##### RANDOM BASELINE #####
setwd("/media/gioca90/University/TiU/Research/Projects/ba_thesis/output/FSCrandom/")
iters = seq(1, 1000, 1)
stats_rnd = NULL
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
  cor.aoa.osc_te_rnd = cor(df.ortho$OSC_te_rnd.z, df.ortho$AoA)
  cor.aoa.osc_ld_rnd = cor(df.ortho$OSC_ld_rnd.z, df.ortho$AoA)

  df.phono$PSC_te_rnd.z = c(scale(box_cox_transf(df.phono$PSC_te_rnd + 0.0001, -10, 10), center = T, scale = T))
  df.phono$PSC_ld_rnd.z = c(scale(box_cox_transf(df.phono$PSC_ld_rnd + 0.0001, -10, 10), center = T, scale = T))
  cor.aoa.psc_te_rnd = cor(df.phono$PSC_te_rnd.z, df.phono$AoA)
  cor.aoa.psc_ld_rnd = cor(df.phono$PSC_ld_rnd.z, df.phono$AoA)

  osc_model.te.rnd <- lm(formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + OSC_te_rnd.z, data = df.ortho)
  deltaAIC.osc_te_rnd = base.ortho.aic - AIC(osc_model.te.rnd)
  
  osc_model.ld.rnd <- lm(formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph + OSC_ld_rnd.z, data = df.ortho)
  deltaAIC.osc_ld_rnd = base.ortho.aic - AIC(osc_model.ld.rnd)
  
  psc_model.te.rnd <- lm(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + PSC_te_rnd.z, data = df.phono)
  deltaAIC.psc_te_rnd = base.phono.aic - AIC(psc_model.te.rnd)
  
  psc_model.ld.rnd <- lm(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph + PSC_ld_rnd.z, data = df.phono)
  deltaAIC.psc_ld_rnd = base.phono.aic - AIC(psc_model.ld.rnd)
  
  df.phono = subset(df.phono, select = -c(PSC_te_rnd, PSC_ld_rnd))
  df.ortho = subset(df.ortho, select = -c(OSC_te_rnd, OSC_ld_rnd))
  
  stats_rnd = rbind(stats_rnd,
                data.frame(i, 
                           cor.aoa.osc_te_rnd, cor.aoa.osc_ld_rnd, cor.aoa.psc_te_rnd, cor.aoa.psc_ld_rnd, 
                           deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd))
}

rm(i, 
   deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd,
   cor.aoa.osc_te_rnd, cor.aoa.osc_ld_rnd, cor.aoa.psc_te_rnd, cor.aoa.psc_ld_rnd,
   psc_model.ld.rnd, psc_model.te.rnd, osc_model.ld.rnd, osc_model.te.rnd, df_rnd)

subtitle_osc_rnd.te = expression(OSC[te-rnd])
subtitle_osc_rnd.ld = expression(OSC[ld-rnd])
subtitle_psc_rnd.te = expression(PSC[te-rnd])
subtitle_psc_rnd.ld = expression(PSC[ld-rnd])

# Correlations
p1.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.osc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.ortho$AoA,df.ortho$OSC_te.z), color = "darkred", size = 2) +
  ggtitle(expression("objective AoA and" ~ OSC[rnd] ~ ": Pearson's r"), subtitle = subtitle_osc_rnd.te)  +
  labs(x = 'r', y = 'freq') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 100)))
(1+sum(stats_rnd$cor.AoA.osc_te_rnd <= cor(df.ortho$AoA,df.ortho$OSC_te.z)))/(length(stats_rnd$cor.AoA.osc_te_rnd)+1)

p2.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.osc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.ortho$AoA,df.ortho$OSC_ld.z), color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_rnd.ld)  +
  labs(x = 'r', y = ' ') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 100)))
(1+sum(stats_rnd$cor.AoA.osc_ld_rnd <= cor(df.ortho$AoA,df.ortho$OSC_ld.z)))/(length(stats_rnd$cor.AoA.osc_ld_rnd)+1)

p3.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.psc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.phono$AoA,df.phono$PSC_te.z), color = "darkred", size = 2) +
  ggtitle(expression("objective AoA and" ~ PSC[rnd] ~ ": Pearson's r"), subtitle = subtitle_psc_rnd.te)  +
  labs(x = 'r', y = 'freq') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 120)))
(1+sum(stats_rnd$cor.AoA.psc_te_rnd <= cor(df.phono$AoA,df.phono$PSC_te.z)))/(length(stats_rnd$cor.AoA.psc_te_rnd)+1)

p4.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.psc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.phono$AoA,df.phono$PSC_ld.z), color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.ld)  +
  labs(x = 'r', y = ' ') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 120)))
(1+sum(stats_rnd$cor.AoA.psc_ld_rnd <= cor(df.phono$AoA,df.phono$PSC_ld.z)))/(length(stats_rnd$cor.AoA.psc_ld_rnd)+1)

grid.arrange(p1.fsc_rnd, p2.fsc_rnd, ncol = 2)
grid.arrange(p3.fsc_rnd, p4.fsc_rnd, ncol = 2)

# DeltaAIC, linear models
range(stats_rnd$deltaAIC.osc_te_rnd)
p1.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.osc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base.ortho.aic - osc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("objective AoA and" ~ OSC[rnd] ~": LM " ~ Delta[AIC]), subtitle = subtitle_osc_rnd.te) +
  labs(x = expression(Delta[AIC]), y = 'freq') +
  coord_cartesian(xlim = (c(-2.5, 20)), ylim = (c(0, 650)))
(1+sum(stats_rnd$deltaAIC.osc_te_rnd >= (base_model.ortho.int.aic - osc_model.te.aic)))/(length(iters)+1)
# 0.000999001

range(stats_rnd$deltaAIC.osc_ld_rnd)
p2.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.osc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base.ortho.aic - osc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_rnd.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')  +
  coord_cartesian(xlim = (c(-2.5, 135)), ylim = (c(0, 650)))
(1+sum(stats_rnd$deltaAIC.osc_ld_rnd >= (base_model.ortho.int.aic - osc_model.ld.aic)))/(length(iters)+1)
# 0.000999001

range(stats_rnd$deltaAIC.psc_te_rnd)
p3.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.psc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base.phono.aic - psc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("objective AoA and" ~ PSC[rnd] ~": LM " ~ Delta[AIC]), subtitle = subtitle_psc_rnd.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq')  +
  coord_cartesian(xlim = (c(-2.5, 18)), ylim = (c(0, 250)))
(1+sum(stats_rnd$deltaAIC.psc_te_rnd >= (base_model.phono.int.aic - psc_model.te.aic)))/(length(iters)+1)
# 0.000999001

range(stats_rnd$deltaAIC.psc_ld_rnd)
p4.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.psc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base.phono.aic - psc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')   +
  coord_cartesian(xlim = (c(-2.5, 260)), ylim = (c(0, 250)))
(1+sum(stats_rnd$deltaAIC.psc_ld_rnd >= (base_model.phono.int.aic - psc_model.ld.aic)))/(length(iters)+1)
# 0.000999001

grid.arrange(p3.fsc_rnd.aic, p4.fsc_rnd.aic, ncol = 2)
grid.arrange(p1.fsc_rnd.aic, p2.fsc_rnd.aic, ncol = 2)


##### VOCABULARY SUBSET BASELINE #####
parent_dir = "/media/gioca90/University/TiU/Research/Projects/ba_thesis/output/FSCsubset/"
setwd(parent_dir)
iters = seq(1, 500, 1)
rates = c(50, 75)
stats_subset = NULL
for (r in rates) {
  print(sprintf("Processing subsets at %s of the whole reference vocabulary", r))
  print(getwd())
  setwd(sprintf('rate%s', r))
  for (i in iters) {
    print(sprintf("Processing iteration %s", i))
    df_subset = read.csv(sprintf("df%s.csv", i), header = T, sep = ',')  
    
    colnames(df_subset)[which(names(df_subset) == "OSC_te")] <- "OSC_te_subset"
    colnames(df_subset)[which(names(df_subset) == "OSC_ld")] <- "OSC_ld_subset"
    colnames(df_subset)[which(names(df_subset) == "PSC_te")] <- "PSC_te_subset"
    colnames(df_subset)[which(names(df_subset) == "PSC_ld")] <- "PSC_ld_subset"
    
    df.phono = merge(df.phono, df_subset[c('word', "PSC_te_subset", "PSC_ld_subset")], by='word')
    df.ortho = merge(df.ortho, df_subset[c('word', "OSC_te_subset", "OSC_ld_subset")], by='word')
    
    df.ortho$OSC_te_subset.z = c(scale(box_cox_transf(df.ortho$OSC_te_subset + 0.0001, -10, 10), center = T, scale = T))
    df.ortho$OSC_ld_subset.z = c(scale(box_cox_transf(df.ortho$OSC_ld_subset + 0.0001, -10, 10), center = T, scale = T))
    
    df.phono$PSC_te_subset.z = c(scale(box_cox_transf(df.phono$PSC_te_subset + 0.0001, -10, 10), center = T, scale = T))
    df.phono$PSC_ld_subset.z = c(scale(box_cox_transf(df.phono$PSC_ld_subset + 0.0001, -10, 10), center = T, scale = T))
    
    osc_model.te.subset <- lm(formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_te_subset.z, data = df.ortho)
    deltaAIC.osc_te_subset = base_model.ortho.int.aic - AIC(osc_model.te.subset)
    
    osc_model.ld.subset <- lm(formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_ld_subset.z, data = df.ortho)
    deltaAIC.osc_ld_subset = base_model.ortho.int.aic - AIC(osc_model.ld.subset)
    
    psc_model.te.subset <- lm(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_te_subset.z, data = df.phono)
    deltaAIC.psc_te_subset = base_model.phono.int.aic - AIC(psc_model.te.subset)
    
    psc_model.ld.subset <- lm(formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_ld_subset.z, data = df.phono)
    deltaAIC.psc_ld_subset = base_model.phono.int.aic - AIC(psc_model.ld.subset)
    
    df.phono = subset(df.phono, select = -c(PSC_te_subset, PSC_ld_subset))
    df.ortho = subset(df.ortho, select = -c(OSC_te_subset, OSC_ld_subset))
    
    stats_subset = rbind(stats_subset,
                  data.frame(i, r,
                             deltaAIC.osc_te_subset, deltaAIC.osc_ld_subset, deltaAIC.psc_te_subset, deltaAIC.psc_ld_subset))
  }
  setwd(parent_dir)
  print(getwd())
}

rm(i, r, 
   deltaAIC.osc_te_subset, deltaAIC.osc_ld_subset, deltaAIC.psc_te_subset, deltaAIC.psc_ld_subset,
   psc_model.ld.subset, psc_model.te.subset, osc_model.ld.subset, osc_model.te.subset)

stats_subset$r <- mapvalues(stats_subset$r, from = c(50, 75), to = c("sampling rate: 50%", "sampling rate: 75%"))

subtitle_osc_subset.te = expression(OSC[te])
subtitle_osc_subset.ld = expression(OSC[ld])
subtitle_psc_subset.te = expression(PSC[te])
subtitle_psc_subset.ld = expression(PSC[ld])

# DeltaAIC, linear models
range(stats_subset$deltaAIC.osc_te_subset)
p1.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.osc_te_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and OSC (subset): LM " ~ Delta[AIC]), subtitle = subtitle_osc_subset.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq') +
  coord_cartesian(xlim = (c(0, 90)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.osc_te_subset >= (base_model.ortho.int.aic - osc_model.te.aic)))/(length(iters)+1)

range(stats_subset$deltaAIC.osc_ld_subset)
p2.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.osc_ld_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_subset.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')  +
  coord_cartesian(xlim = (c(0, 215)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.osc_ld_subset >= (base_model.ortho.int.aic - osc_model.ld.aic)))/(length(iters)+1)

range(stats_subset$deltaAIC.psc_te_subset)
p3.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.psc_te_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - psc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("objective AoA and PSC (subset): LM " ~ Delta[AIC]), subtitle = subtitle_psc_subset.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq')  +
  coord_cartesian(xlim = (c(0, 90)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
stats_subset_75 = stats_subset[stats_subset$r == 'sampling rate: 75%',]
stats_subset_50 = stats_subset[stats_subset$r == 'sampling rate: 50%',]
(1+sum(stats_subset_50$deltaAIC.psc_te_subset >= (base_model.phono.int.aic - psc_model.te.aic)))/(500+1)
(1+sum(stats_subset_75$deltaAIC.psc_te_subset >= (base_model.phono.int.aic - psc_model.te.aic)))/(500+1)

range(stats_subset$deltaAIC.psc_ld_subset)
p4.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.psc_ld_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base.phono.aic - psc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_subset.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')   +
  coord_cartesian(xlim = (c(0, 325)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
(1+sum(stats_subset_50$deltaAIC.psc_ld_subset >= (base_model.phono.int.aic - psc_model.ld.aic)))/(500+1)
(1+sum(stats_subset_75$deltaAIC.psc_ld_subset >= (base_model.phono.int.aic - psc_model.ld.aic)))/(500+1)

grid.arrange(p3.fsc_subset.aic, p4.fsc_subset.aic, ncol = 2)
grid.arrange(p1.fsc_subset.aic, p2.fsc_subset.aic, ncol = 2)

range(stats$deltaAIC.psc_ld_subset)


##### PCA FOR NEIGHBORHOOD MEASURES - ORTHOGRAPHY #####
neighborhoods_ld.ortho.pca <- prcomp(df.ortho[c("OSC_ld.z", "old.z", "snd.z")], center = TRUE, scale. = TRUE)

plot(cumsum(neighborhoods_ld.ortho.pca$sdev^2 / sum(neighborhoods_ld.ortho.pca$sdev^2)), type="b", col='steelblue')

pc_ortho = df.ortho[c("OSC_ld.z", "old.z", "snd.z")]
colnames(pc_ortho) <- c("OSCld", "OLD20", "SND")
corrplot(cor(neighborhoods_ld.ortho.pca$x[,1:3], pc_ortho), 
         title = "Correlation matrix: PCs",
         mar = c(2,2,2,2),
         method = 'number',
         tl.col = 'black',
         tl.srt = 30)
rm(pc_ortho)

df.ortho['PC1_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,1]
df.ortho['PC2_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,2]
df.ortho['PC3_ld_ortho'] = neighborhoods_ld.ortho.pca$x[,3]


pca_ld_model.ortho <- lm(
  formula = AoA ~ n_chars.z + concr.z + val.z + freq.z + morph + PC1_ld_ortho + PC2_ld_ortho + PC3_ld_ortho, data = df.ortho
)
pca_ld_model.ortho.aic = AIC(pca_ld_model.ortho)
# 38302.7
summary(pca_ld_model.ortho)$r.sq
# 0.4440944
summary(pca_ld_model.ortho)
# all principal components are significant, but PC2 and PC3 only slightly
# PC1 high when OSC high, OLD low, SND low: negative effect on AoA
# PC2 high when OLD higher, SND low: positive effect on AoA
# PC3 hgh when OSC higher, OLD higher, SND slightly higher: negative effect on AoA
base.ortho.aic - pca_ld_model.ortho.aic


##### PCA FOR NEIGHBORHOOD MEASURES - PHONOLOGY #####
neighborhoods_ld.phono.pca <- prcomp(df.phono[c("PSC_ld.z", "pnd.z", "snd.z")], center = TRUE, scale. = TRUE)

plot(cumsum(neighborhoods_ld.phono.pca$sdev^2 / sum(neighborhoods_ld.phono.pca$sdev^2)), type="b", col='steelblue')

pc_phono = df.phono[c("PSC_ld.z", "pnd.z", "snd.z")]
colnames(pc_phono) <- c("PSCld", "PND", "SND")
corrplot(cor(neighborhoods_ld.phono.pca$x[,1:3], pc_phono), 
         title = "Correlation matrix: PCs",
         mar = c(2,2,2,2),
         method = 'number',
         tl.col = 'black',
         tl.srt = 30)
rm(pc_phono)

df.phono['PC1_ld_phono'] = neighborhoods_ld.phono.pca$x[,1]
df.phono['PC2_ld_phono'] = neighborhoods_ld.phono.pca$x[,2]
df.phono['PC3_ld_phono'] = neighborhoods_ld.phono.pca$x[,3]


pca_ld_model.phono <- lm(
  formula = AoA ~ n_phons.z + concr.z + val.z + freq.z + morph + PC1_ld_phono + PC2_ld_phono + PC3_ld_phono, data = df.phono
)
pca_ld_model.phono.aic = AIC(pca_ld_model.phono)
# 30194.1
summary(pca_ld_model.phono)$r.sq
# 0.4788061
summary(pca_ld_model.phono)
# only PC1 and PC3 are significant, PC2 < 0.1
# PC1 high when PSCld high, PND high, SND lower: negative effect on AoA
# PC3 high when PSCld low, PND higher, SND slightly lower: positive effect on AoA
base.phono.aic - pca_ld_model.phono.aic
# 248.3814

