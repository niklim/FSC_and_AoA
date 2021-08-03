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
library(psych)
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

#transform morph into factor
data$morph <- factor(data$morph,    
                     level=c(0,1),    
                     labels=c("Mono","Poly"))


#get overview of the variables
head(data, 5)
summary(data)
describe(data)


### Transform variables into normal shape and comparable units (SD) ###
df.ortho <- data[,FALSE]
df.phono <- data[,FALSE]

df.ortho$word = data$word
df.phono$word = data$word

df.ortho$morph = data$morph
df.phono$morph = data$morph

df.ortho$aoa = data$aoa
df.phono$aoa = data$aoa

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

df.phono$n_phons.z = c(scale(box_cox_transf(data$n_phon, -10, 10), center = T, scale = T))
df.phono$PSC_te.z = c(scale(box_cox_transf(data$PSC_te + 0.0001, -10, 10), center = T, scale = T))
df.phono$PSC_ld.z = c(scale(box_cox_transf(data$PSC_ld, -10, 10), center = T, scale = T))
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

rm(df.mald, data, iconicity)

##### DESCRIPTIVE PLOTS ####
#get histograms of variables
hist(data$aoa)

hist(data$length)
hist(df$n_chars.z)

hist(data$OSC_te, breaks = 30, main = "Histogram of target embedded OSC Values", xlab = "OSC Value")
hist(df$OSC_te.z, breaks = 30, main = "Histogram of z-transformed target embedded OSC Values", xlab = "OSC Values")

hist(data$OSC_ld, breaks = 30)
hist(df$OSC_ld.z, breaks = 30)

hist(data$PSC_te, breaks = 30, main = "Histogram of target embedded PSC Values", xlab = "OSC Value")
hist(df$PSC_te.z, breaks = 30, main = "Histogram of z-transformed target embedded PSC Values", xlab = "OSC Value")

hist(data$PSC_ld, breaks = 30)
hist(df$PSC_ld.z, breaks = 30)

hist(data$concr)
hist(df$concr.z)

hist(data$val)
hist(df$val.z)

hist(data$freq)
hist(df$freq.z)

hist(data$snd)
hist(df$snd.z)

mono <- filter(df, morph == "Mono")
poly <- filter(df, morph == "Poly")

hist(mono$OSC_te.z, breaks = 30, main = "target embedded OSC Values for Monomorphemic Words", xlab = "OSC Values" )
hist(poly$OSC_te.z, breaks = 30, main = "target embedded OSC Values for Polymorphemic Words", xlab = "OSC Values")

hist(mono$PSC_te.z, breaks = 30, main = "target embedded PSC Values for Monomorphemic Words", xlab = "PSC Values" )
hist(poly$PSC_te.z, breaks = 30, main = "target embedded PSC Values for Polymorphemic Words", xlab = "PSC Values")

hist(mono$OSC_ld.z, breaks = 30, main = "Levenshtein OSC Values for Monomorphemic Words", xlab = "OSC Values" )
hist(poly$OSC_ld.z, breaks = 30, main = "Levenshtein OSC Values for Polymorphemic Words", xlab = "OSC Values")

hist(mono$PSC_ld.z, breaks = 30, main = "Levenshtein PSC Values for Monomorphemic Words", xlab = "PSC Values" )
hist(poly$PSC_ld.z, breaks = 30, main = "Levenshtein PSC Values for Polymorphemic Words", xlab = "PSC Values")

##### CORRELATION MATRIX #####
cor.test(df.phono.icon$aoa, df.phono.icon$iconicity.z)
# -0.07242057 [-0.11829096, -0.02624176]
# t = -3.0746, df = 1793, p-value = 0.002139

cor.test(df.ortho.icon$aoa, df.ortho.icon$iconicity.z)
# -0.04624576 [-0.089550043, -0.002766965]
# t = -2.0859, df = 2030, p-value = 0.03712

cor.test(df.phono.icon$PSC_ld.z, df.phono.icon$iconicity.z)
# 0.1633664 [0.1179915, 0.2080605]
# t = 7.0118, df = 1793, p-value = 3.321e-12

cor.test(df.phono.icon$PSC_te.z, df.phono.icon$iconicity.z)
# -0.008196443 [-0.05444254  0.03808474]
# t = -0.34708, df = 1793, p-value = 0.7286

cor.test(df.ortho.icon$OSC_ld.z, df.ortho.icon$iconicity.z)
# 0.2105296 [0.1685886 0.2517096]
# t = 9.703, df = 2030, p-value < 2.2e-16

cor.test(df.ortho.icon$OSC_te.z, df.ortho.icon$iconicity.z)
# -0.06371315 [-0.10690134 -0.02028499]
# t = -2.8765, df = 2030, p-value = 0.004063

phono <- df.phono[c("aoa", "PSC_te.z", "PSC_ld.z", "concr.z", "val.z", "n_phons.z", "freq.z", "pnd.z", "snd.z")]
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
             main = 'Phonological variables - Subjective AoA')
rm(phono)

ortho <- df.ortho[c("aoa", "OSC_te.z", "OSC_ld.z", "concr.z", "val.z", "n_chars.z", "freq.z", "old.z", "snd.z")]
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
             main = "Orthographic variables - Subjective AoA")
rm(ortho)


##### MULTIPLE REGRESSION - ORTHOGRAPHY #####
# Base Model
base_model.ortho <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z + snd.z + morph, data = df.ortho
  )
summary(base_model.ortho)
VIF(base_model.ortho)
base.ortho.aic <- AIC(base_model.ortho)

# Base Model with snd*old interaction
base_model.ortho.int <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z * snd.z, data = df.ortho
  )
base_model.ortho.int.aic = AIC(base_model.ortho.int)
summary(base_model.ortho.int)$r.sq
base.ortho.aic - base_model.ortho.int.aic
# 10.05795

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
quantile(df.ortho$old.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
base_model.ortho.int.df = as.data.frame(
  ggpredict(
    base_model.ortho.int, terms = c(
      "old.z[-1.872, -1.506, -0.724, -0.482, -0.26, -0.122, 0.310, 0.632, 0.872, 1.335, 3.477]",
      "snd.z[-1.61132041, -0.02946804, 1.65420042]")
    )
  )
base_model.ortho.int.df$group <- mapvalues(
  base_model.ortho.int.df$group, 
  from = c("-1.61132041", "-0.02946804", "1.65420042"), 
  to = c("5%", "50%", "95%")
  )

ggplot(base_model.ortho.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('OLD') +
  ylab('AoA') +
  ggtitle(expression("Interaction of OLD20 and SND on subjective AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5, 3)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


# Iconicity control
icon_model.ortho <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + iconicity.z, data = df.ortho.icon
  )
icon_model.ortho.aic = AIC(icon_model.ortho)
summary(icon_model.ortho)$r.sq
summary(icon_model.ortho)
VIF(icon_model.ortho)

# target embedded neighbors
osc_model.te <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_te.z, data = df.ortho
  )
osc_model.te.aic = AIC(osc_model.te)
summary(osc_model.te)$r.sq
summary(osc_model.te)
VIF(osc_model.te)
base_model.ortho.int.aic - osc_model.te.aic
# 5.898347

# target embedded neighbors over iconicity
osc_icon_model.te <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + iconicity.z + OSC_te.z, data = df.ortho.icon
  )
osc_icon_model.te.aic = AIC(osc_icon_model.te)
summary(osc_icon_model.te)$r.sq
summary(osc_icon_model.te)
VIF(osc_icon_model.te)
icon_model.ortho.aic - osc_icon_model.te.aic
# 7.508975

# OSCte*morph interaction
osc_by_morph_model.te <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph*OSC_te.z, data = df.ortho
  )
osc_by_morph_model.te.aic = AIC(osc_by_morph_model.te)
summary(osc_by_morph_model.te)$r.sq
summary(osc_by_morph_model.te)
osc_model.te.aic - osc_by_morph_model.te.aic
# 80.80363

# Levenshtein neighbors
osc_model.ld <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_ld.z, data = df.ortho
  )
osc_model.ld.aic = AIC(osc_model.ld)
summary(osc_model.ld)$r.sq
summary(osc_model.ld)
base_model.ortho.int.aic - osc_model.ld.aic
# 266.8096

# Levenshtein neighbors over iconicity
osc_icon_model.ld <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + iconicity.z + OSC_ld.z, data = df.ortho.icon
  )
osc_icon_model.ld.aic = AIC(osc_icon_model.ld)
summary(osc_icon_model.ld)$r.sq
summary(osc_icon_model.ld)
VIF(osc_icon_model.ld)
icon_model.ortho.aic - osc_icon_model.ld.aic
# 84.81813

# OSCld*morph interaction
osc_by_morph_model.ld <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph*OSC_ld.z, data = df.ortho
  )
osc_by_morph_model.ld.aic = AIC(osc_by_morph_model.ld)
summary(osc_by_morph_model.ld)$r.sq
summary(osc_by_morph_model.ld)
osc_model.ld.aic - osc_by_morph_model.ld.aic
# 5.240545


##### MULTIPLE REGRESSION - PHONOLOGY #####
# Base Model
base_model.phono <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + snd.z + morph, data = df.phono
  )
summary(base_model.phono)$r.sq
summary(base_model.phono)
VIF(base_model.phono)
base.phono.aic <- AIC(base_model.phono)

# Base Model with pnd*snd interaction
base_model.phono.int <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph, data = df.phono
  )
summary(base_model.phono.int)$r.sq
summary(base_model.phono.int)
VIF(base_model.phono.int)
base_model.phono.int.aic <- AIC(base_model.phono.int)
base.phono.aic - base_model.phono.int.aic

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
quantile(df.phono$pnd.z, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
base_model.phono.int.df = as.data.frame(
  ggpredict(
    base_model.phono.int, terms = c(
      "pnd.z[-2.126, -0.820, -0.820, -0.59, -0.325, -0.082, 0.295, 0.605, 0.935, 1.317, 2.196]",
      "snd.z[-1.63825643, -0.03833438, 1.61215118]")
  )
)
base_model.phono.int.df$group <- mapvalues(
  base_model.phono.int.df$group, 
  from = c("-1.63825643", "-0.03833438", "1.61215118"), 
  to = c("5%", "50%", "95%")
)

ggplot(base_model.phono.int.df, aes(x = x, y = predicted, colour = group, linetype = group)) +  
  geom_line(alpha = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab('PND') +
  ylab('AoA') +
  ggtitle(expression("Interaction of PND and SND on subjective AoA")) +
  labs(colour = "SND", fill = "SND", linetype = 'SND') +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# Iconicity control
icon_model.phono <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + iconicity.z, data = df.phono.icon
  )
icon_model.phono.aic = AIC(icon_model.phono)
summary(icon_model.phono)$r.sq
summary(icon_model.phono)
VIF(icon_model.phono)

# target embedded neighbors
psc_model.te <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_te.z, data = df.phono
  )
psc_model.te.aic = AIC(psc_model.te)
summary(psc_model.te)$r.sq
summary(psc_model.te)
VIF(psc_model.te)
base_model.phono.int.aic - psc_model.te.aic
# 4.867679

# target embedded neighbors over iconicity
psc_icon_model.te <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + iconicity.z + PSC_te.z, data = df.phono.icon
  )
psc_icon_model.te.aic = AIC(psc_icon_model.te)
summary(psc_icon_model.te)$r.sq
summary(psc_icon_model.te)
VIF(psc_icon_model.te)
icon_model.phono.aic - psc_icon_model.te.aic

# PSCte*morph interaction
psc_by_morph_model.te <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph*PSC_te.z, data = df.phono
  )
psc_by_morph_model.te.aic = AIC(psc_by_morph_model.te)
summary(psc_by_morph_model.te)$r.sq
summary(psc_by_morph_model.te)
psc_model.te.aic - psc_by_morph_model.te.aic

# Levenshtein neighbors
psc_model.ld <- lm(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_ld.z, data = df.phono)
psc_model.ld.aic = AIC(psc_model.ld)
summary(psc_model.ld)$r.sq
summary(psc_model.ld)
base_model.phono.int.aic - psc_model.ld.aic

# Levenshtein neighbors over iconicity
psc_icon_model.ld <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + iconicity.z + PSC_ld.z, data = df.phono.icon
  )
psc_icon_model.ld.aic = AIC(psc_icon_model.ld)
summary(psc_icon_model.ld)$r.sq
summary(psc_icon_model.ld)
VIF(psc_icon_model.ld)
icon_model.phono.aic - psc_icon_model.ld.aic

# PSCld*morph interaction
psc_by_morph_model.ld <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph*PSC_ld.z, data = df.phono
  )
psc_by_morph_model.ld.aic = AIC(psc_by_morph_model.ld)
summary(psc_by_morph_model.ld)$r.sq
summary(psc_by_morph_model.ld)
psc_model.ld.aic - psc_by_morph_model.ld.aic

##### Interaction FSC*SND #####
### replication iconicity*SND interaction with FSC measures ###
cor.test(df.ortho.icon$iconicity.z, df.ortho.icon$snd.z)
cor.test(df.phono.icon$iconicity.z, df.phono.icon$snd.z)

cor.test(df.ortho.icon$iconicity.z, df.ortho.icon$OSC_ld.z)
cor.test(df.phono.icon$iconicity.z, df.phono.icon$PSC_ld.z)

cor.test(df.ortho$OSC_ld.z, df.ortho$snd.z)
cor.test(df.phono$PSC_ld.z, df.phono$snd.z)

# predict systematicity (FSC-ld)
OSC_ld.lm <- lm(formula = OSC_ld.z ~ n_chars.z + freq.z + old.z + snd.z + aoa, data = df.ortho)
summary(OSC_ld.lm)

PSC_ld.lm <- lm(formula = PSC_ld.z ~ n_phons.z + freq.z + pnd.z + snd.z + aoa, data = df.phono)
summary(PSC_ld.lm)


# 2-way interaction: SND*OSC
osc_by_snd_model.ld <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_ld.z * snd.z, data = df.ortho
  )
osc_by_snd_model.ld.aic = AIC(osc_by_snd_model.ld)
summary(osc_by_snd_model.ld)$r.sq
summary(osc_by_snd_model.ld)
osc_model.ld.aic - osc_by_snd_model.ld.aic

quantile(df.ortho$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.ortho = data.frame(
  ggpredict(osc_by_snd_model.ld, terms = c("OSC_ld.z","snd.z[-1.61132041, -0.02946804, 1.65420042]"))
)
effects.2way_int.ortho$group <- mapvalues(
  effects.2way_int.ortho$group, 
  from = c("-1.61132041", "-0.02946804", "1.65420042"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.ortho, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('OSC'[ld])) +
  ylab('AoA') +
  ggtitle(expression("2way interaction of SND and" ~ OSC[ld] ~ "on AoA")) +
  labs(colour = "SND", fill = "SND")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# 2-way interaction: SND*PSC
psc_by_snd_model.ld <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_ld.z * snd.z, data = df.phono
)
psc_by_snd_model.ld.aic = AIC(psc_by_snd_model.ld)
summary(psc_by_snd_model.ld)$r.sq
summary(psc_by_snd_model.ld)
psc_model.ld.aic - psc_by_snd_model.ld.aic

quantile(df.phono$snd.z, c(0.05, 0.5, 0.95))
effects.2way_int.phono = data.frame(
  ggpredict(psc_by_snd_model.ld, terms = c("PSC_ld.z","snd.z[-1.63825643, -0.03833438, 1.61215118]"))
)
effects.2way_int.phono$group <- mapvalues(
  effects.2way_int.phono$group, 
  from = c("-1.63825643", "-0.03833438", "1.61215118"), 
  to = c("5%", "50%", "95%")
)

ggplot(effects.2way_int.phono, aes(x = x, y = predicted, colour = group)) +  
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  xlab(expression('PSC'[ld])) +
  ylab('AoA') +
  ggtitle(expression("2way interaction of SND and" ~ PSC[ld] ~ "on AoA")) +
  labs(colour = "SND", fill = "SND")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


##### Random Forest Regression Models #####

### Random Forest Base ###
rf_base.ortho <- ranger(formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z + snd.z, data = df.ortho, mtry = 5/3, importance = "impurity" )
print(rf_base.ortho)
vip(rf_base.ortho, aesthetics = list(col = "darkblue", fill = "darkblue"))
rf_base.ortho.rsq = rf_base.ortho$r.squared

rf_base.phono <- ranger(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + morph + pnd.z + snd.z, data = df.phono, mtry = 5/3, importance = "impurity" )
print(rf_base.phono)
vip(rf_base.phono, aesthetics = list(col = "darkblue", fill = "darkblue"))
rf_base.phono.rsq = rf_base.phono$r.squared


### Random Forest OSCte ###
rf_OSC_te <- ranger(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + OSC_te.z, data = df.ortho, mtry = 6/3, importance = "impurity"
  )
print(rf_OSC_te)
names(rf_OSC_te$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "OSC[te]"
  ))
rf_OSC_te.vip_df = data.frame(
  var=names(rf_OSC_te$variable.importance), Importance=rf_OSC_te$variable.importance, row.names=NULL
)
subtitle_osc.te = parse(text = "OSC[te]")
p1.rf = ggplot(data = rf_OSC_te.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs) - subjective AoA", subtitle = subtitle_osc.te) +
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
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + OSC_ld.z, data = df.ortho, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_ld)
names(rf_OSC_ld$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "OSC[ld]"
))
rf_OSC_ld.vip_df = data.frame(
  var=names(rf_OSC_ld$variable.importance), Importance=rf_OSC_ld$variable.importance, row.names=NULL
)
subtitle_osc.ld = parse(text = "OSC[ld]")
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
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + PSC_te.z, data = df.phono, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_te)
names(rf_PSC_te$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "PSC[te]"
))
rf_PSC_te.vip_df = data.frame(
  var=names(rf_PSC_te$variable.importance), Importance=rf_PSC_te$variable.importance, row.names=NULL
)
subtitle_psc.te = parse(text = "PSC[te]")
p3.rf = ggplot(data = rf_PSC_te.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs) - subjective AoA", subtitle = subtitle_psc.te) +
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
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + PSC_ld.z, data = df.phono, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_ld)
names(rf_PSC_ld$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "PSC[ld]"
))
rf_PSC_ld.vip_df = data.frame(
  var=names(rf_PSC_ld$variable.importance), Importance=rf_PSC_ld$variable.importance, row.names=NULL
)
subtitle_psc.ld = parse(text = "PSC[ld]")
p4.rf = ggplot(data = rf_PSC_ld.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle(" ", subtitle = subtitle_psc.ld) +
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

grid.arrange(p3.rf, p4.rf, ncol = 1)
grid.arrange(p1.rf, p2.rf, ncol = 1)



### Random Forest Iconicity Base ###
rf_base.ortho.icon <- ranger(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + morph + old.z + snd.z + iconicity.z, data = df.ortho.icon, mtry = 5/3, importance = "impurity"
  )
print(rf_base.ortho.icon)
rf_base.ortho.icon.rsq = rf_base.ortho.icon$r.squared

rf_base.phono.icon <- ranger(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + morph + pnd.z + snd.z + iconicity.z, data = df.phono.icon, mtry = 5/3, importance = "impurity"
  )
print(rf_base.phono.icon)
rf_base.phono.icon.rsq = rf_base.phono.icon$r.squared


### Random Forest Iconicity + OSCte ###
rf_OSC_te.icon <- ranger(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + iconicity.z + OSC_te.z, 
  data = df.ortho.icon, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_te.icon)
names(rf_OSC_te.icon$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "Icon", "OSC[te]"
))
rf_OSC_te.icon.vip_df = data.frame(
  var=names(rf_OSC_te.icon$variable.importance), Importance=rf_OSC_te.icon$variable.importance, row.names=NULL
)
subtitle_osc.te.icon = parse(text = "icon + OSC[te]")
p1.rf.icon = ggplot(data = rf_OSC_te.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs w/ iconicity) - subjective AoA", subtitle = subtitle_osc.te.icon) +
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
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z + morph + snd.z + iconicity.z + OSC_ld.z, 
  data = df.ortho.icon, mtry = 6/3, importance = "impurity"
)
print(rf_OSC_ld.icon)
names(rf_OSC_ld.icon$variable.importance) <- parse(text = c(
  "Len[ortho]", "Concr",  "Val", "Freq", "OLD20", "Morph", "SND", "Icon", "OSC[ld]"
))
rf_OSC_ld.icon.vip_df = data.frame(
  var=names(rf_OSC_ld.icon$variable.importance), Importance=rf_OSC_ld.icon$variable.importance, row.names=NULL
)
subtitle_osc.ld.icon = parse(text = "icon + OSC[ld]")
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
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + iconicity.z + PSC_te.z, 
  data = df.phono.icon, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_te.icon)
names(rf_PSC_te.icon$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "Icon", "PSC[te]"
))
rf_PSC_te.icon.vip_df = data.frame(
  var=names(rf_PSC_te.icon$variable.importance), Importance=rf_PSC_te.icon$variable.importance, row.names=NULL
)
subtitle_psc.te.icon = parse(text = "icon + PSC[te]")
p3.rf.icon = ggplot(data = rf_PSC_te.icon.vip_df, aes(y = Importance, x = reorder(var, Importance))) +
  geom_col(color = "steelblue", fill = "steelblue") +
  ggtitle("Variable Importance (RFs w/ iconicity) - subjective AoA", subtitle = subtitle_psc.te.icon) +
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
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z + morph + snd.z + iconicity.z + PSC_ld.z, 
  data = df.phono.icon, mtry = 6/3, importance = "impurity"
)
print(rf_PSC_ld.icon)
names(rf_PSC_ld.icon$variable.importance) <- parse(text = c(
  "Len[phono]", "Concr",  "Val", "Freq", "PND", "Morph", "SND", "Icon", "PSC[ld]"
))
rf_PSC_ld.icon.vip_df = data.frame(
  var=names(rf_PSC_ld.icon$variable.importance), Importance=rf_PSC_ld.icon$variable.importance, row.names=NULL
)
subtitle_psc.ld.icon = parse(text = "icon + PSC[ld]")
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

grid.arrange(p3.rf.icon, p4.rf.icon, ncol = 1)
grid.arrange(p1.rf.icon, p2.rf.icon, ncol = 1)


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
  cor.aoa.osc_te_rnd = cor(df.ortho$OSC_te_rnd.z, df.ortho$aoa)
  cor.aoa.osc_ld_rnd = cor(df.ortho$OSC_ld_rnd.z, df.ortho$aoa)
  
  df.phono$PSC_te_rnd.z = c(scale(box_cox_transf(df.phono$PSC_te_rnd + 0.0001, -10, 10), center = T, scale = T))
  df.phono$PSC_ld_rnd.z = c(scale(box_cox_transf(df.phono$PSC_ld_rnd + 0.0001, -10, 10), center = T, scale = T))
  cor.aoa.psc_te_rnd = cor(df.phono$PSC_te_rnd.z, df.phono$aoa)
  cor.aoa.psc_ld_rnd = cor(df.phono$PSC_ld_rnd.z, df.phono$aoa)
  
  osc_model.te.rnd <- lm(formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_te_rnd.z, data = df.ortho)
  deltaAIC.osc_te_rnd = base_model.ortho.int.aic - AIC(osc_model.te.rnd)
  
  osc_model.ld.rnd <- lm(formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_ld_rnd.z, data = df.ortho)
  deltaAIC.osc_ld_rnd = base_model.ortho.int.aic - AIC(osc_model.ld.rnd)
  
  psc_model.te.rnd <- lm(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_te_rnd.z, data = df.phono)
  deltaAIC.psc_te_rnd = base_model.phono.int.aic - AIC(psc_model.te.rnd)
  
  psc_model.ld.rnd <- lm(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_ld_rnd.z, data = df.phono)
  deltaAIC.psc_ld_rnd = base_model.phono.int.aic - AIC(psc_model.ld.rnd)
  
  df.phono = subset(df.phono, select = -c(PSC_te_rnd, PSC_ld_rnd))
  df.ortho = subset(df.ortho, select = -c(OSC_te_rnd, OSC_ld_rnd))

  stats_rnd = rbind(stats_rnd,
                data.frame(i, 
                           cor.aoa.osc_te_rnd, cor.aoa.osc_ld_rnd, cor.aoa.psc_te_rnd, cor.aoa.psc_ld_rnd, 
                           deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd))
}

rm(i, cor.aoa.osc_te_rnd, cor.aoa.osc_ld_rnd, cor.aoa.psc_te_rnd, cor.aoa.psc_ld_rnd, 
   deltaAIC.osc_te_rnd, deltaAIC.osc_ld_rnd, deltaAIC.psc_te_rnd, deltaAIC.psc_ld_rnd)

subtitle_osc_rnd.te = expression(OSC[te-rnd])
subtitle_osc_rnd.ld = expression(OSC[ld-rnd])
subtitle_psc_rnd.te = expression(PSC[te-rnd])
subtitle_psc_rnd.ld = expression(PSC[ld-rnd])

# Correlations
p1.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.osc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.ortho$aoa,df.ortho$OSC_te.z), color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and" ~ OSC[rnd] ~ ": Pearson's r"), subtitle = subtitle_osc_rnd.te)  +
  labs(x = 'r', y = 'freq') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 100)))
(1+sum(stats_rnd$cor.aoa.osc_te_rnd <= cor(df.ortho$aoa,df.ortho$OSC_te.z)))/(length(stats_rnd$cor.aoa.osc_te_rnd)+1)

p2.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.osc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.ortho$aoa,df.ortho$OSC_ld.z), color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_rnd.ld)  +
  labs(x = 'r', y = ' ') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 100)))
(1+sum(stats_rnd$cor.aoa.osc_ld_rnd <= cor(df.ortho$aoa,df.ortho$OSC_ld.z)))/(length(stats_rnd$cor.aoa.osc_ld_rnd)+1)

p3.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.psc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.phono$aoa,df.phono$PSC_te.z), color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and" ~ PSC[rnd] ~ ": Pearson's r"), subtitle = subtitle_psc_rnd.te)  +
  labs(x = 'r', y = 'freq') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 130)))
(1+sum(stats_rnd$cor.aoa.psc_te_rnd <= cor(df.phono$aoa,df.phono$PSC_te.z)))/(length(stats_rnd$cor.aoa.psc_te_rnd)+1)

p4.fsc_rnd = ggplot(data = stats_rnd, aes(x=cor.aoa.psc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = cor(df.phono$aoa,df.phono$PSC_ld.z), color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.ld)  +
  labs(x = 'r', y = ' ') +
  coord_cartesian(xlim = (c(-0.55, 0.1)), ylim = (c(0, 130)))
(1+sum(stats_rnd$cor.aoa.psc_ld_rnd <= cor(df.phono$aoa,df.phono$PSC_ld.z)))/(length(stats_rnd$cor.aoa.psc_ld_rnd)+1)

grid.arrange(p1.fsc_rnd, p2.fsc_rnd, ncol = 2)
grid.arrange(p3.fsc_rnd, p4.fsc_rnd, ncol = 2)


# DeltaAIC, linear models
p1.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.osc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and" ~ OSC[rnd] ~ ": LM " ~ Delta[AIC]), subtitle = subtitle_osc_rnd.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq') +
  coord_cartesian(xlim = (c(-2.5, 15)), ylim = (c(0, 550)))
(1+sum(stats_rnd$deltaAIC.osc_te_rnd >= (base_model.ortho.int.aic - osc_model.te.aic)))/(length(iters)+1)

p2.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.osc_ld_rnd)) +
    geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
    geom_vline(xintercept = base_model.ortho.int.aic - osc_model.ld.aic, color = "darkred", size = 2) +
    ggtitle(" ", subtitle = subtitle_osc_rnd.ld)  +
    labs(x = expression(Delta[AIC]), y = ' ')  +
  coord_cartesian(xlim = (c(-2.5, 450)), ylim = (c(0, 550)))
(1+sum(stats_rnd$deltaAIC.osc_ld_rnd >= (base_model.ortho.int.aic - osc_model.ld.aic)))/(length(iters)+1)

p3.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.psc_te_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - psc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and" ~ PSC[rnd] ~ ": LM " ~ Delta[AIC]), subtitle = subtitle_psc_rnd.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq')  +
  coord_cartesian(xlim = (c(-2.5, 15)), ylim = (c(0, 350)))
(1+sum(stats_rnd$deltaAIC.psc_te_rnd >= (base_model.phono.int.aic - psc_model.te.aic)))/(length(iters)+1)

p4.fsc_rnd.aic = ggplot(data = stats_rnd, aes(x=deltaAIC.psc_ld_rnd)) +
  geom_histogram(bins = 100, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - psc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_rnd.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ')   +
  coord_cartesian(xlim = (c(-2.5, 450)), ylim = (c(0, 350)))
(1+sum(stats_rnd$deltaAIC.psc_ld_rnd >= (base_model.phono.int.aic - psc_model.ld.aic)))/(length(iters)+1)

grid.arrange(p1.fsc_rnd.aic, p2.fsc_rnd.aic, ncol = 2)
grid.arrange(p3.fsc_rnd.aic, p4.fsc_rnd.aic, ncol = 2)

range(stats_rnd$deltaAIC.psc_ld_rnd)


##### EFFECT OF REFERENCE VOCABULARY #####
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
    
    osc_model.te.subset <- lm(formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_te_subset.z, data = df.ortho)
    deltaAIC.osc_te_subset = base_model.ortho.int.aic - AIC(osc_model.te.subset)
    
    osc_model.ld.subset <- lm(formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + old.z*snd.z + morph + OSC_ld_subset.z, data = df.ortho)
    deltaAIC.osc_ld_subset = base_model.ortho.int.aic - AIC(osc_model.ld.subset)
    
    psc_model.te.subset <- lm(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_te_subset.z, data = df.phono)
    deltaAIC.psc_te_subset = base_model.phono.int.aic - AIC(psc_model.te.subset)
    
    psc_model.ld.subset <- lm(formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + pnd.z*snd.z + morph + PSC_ld_subset.z, data = df.phono)
    deltaAIC.psc_ld_subset = base_model.phono.int.aic - AIC(psc_model.ld.subset)
    
    df.phono = subset(df.phono, select = -c(PSC_te_subset, PSC_ld_subset))
    df.ortho = subset(df.ortho, select = -c(OSC_te_subset, OSC_ld_subset))
    
    stats_subset = rbind(stats_subset,
                  data.frame(i, r, deltaAIC.osc_te_subset, deltaAIC.osc_ld_subset, deltaAIC.psc_te_subset, deltaAIC.psc_ld_subset))
  }
  setwd(parent_dir)
  print(getwd())
}

rm(i, r, 
   deltaAIC.osc_te_subset, deltaAIC.osc_ld_subset, deltaAIC.psc_te_subset, deltaAIC.psc_ld_subset)

subtitle_osc_subset.te = expression(OSC[te])
subtitle_osc_subset.ld = expression(OSC[ld])
subtitle_psc_subset.te = expression(PSC[te])
subtitle_psc_subset.ld = expression(PSC[ld])

stats_subset$r <- mapvalues(stats_subset$r, from = c(50, 75), to = c("sampling rate: 50%", "sampling rate: 75%"))

# DeltaAIC, linear models
range(stats_subset$deltaAIC.osc_te_subset)
p1.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.osc_te_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and OSC (subset): LM " ~ Delta[AIC]), subtitle = subtitle_osc_subset.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq') +
  coord_cartesian(xlim = (c(0, 60)), ylim = (c(0, 50))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.osc_te_subset >= (base_model.ortho.int.aic - osc_model.te.aic)))/(length(iters)+1)

range(stats_subset$deltaAIC.osc_ld_subset)
p2.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.osc_ld_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.ortho.int.aic - osc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_osc_subset.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ') +
  coord_cartesian(xlim = (c(0, 400)), ylim = (c(0, 50))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.osc_ld_subset >= (base_model.ortho.int.aic - osc_model.ld.aic)))/(length(iters)+1)

range(stats_subset$deltaAIC.psc_te_subset)
p3.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.psc_te_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - psc_model.te.aic, color = "darkred", size = 2) +
  ggtitle(expression("subjective AoA and PSC (subset): LM " ~ Delta[AIC]), subtitle = subtitle_psc_subset.te)  +
  labs(x = expression(Delta[AIC]), y = 'freq')  +
  coord_cartesian(xlim = (c(0, 60)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.psc_te_subset >= (base_model.phono.int.aic - psc_model.te.aic)))/(length(iters)+1)

range(stats_subset$deltaAIC.psc_ld_subset)
p4.fsc_subset.aic = ggplot(data = stats_subset, aes(x=deltaAIC.psc_ld_subset)) +
  geom_histogram(bins = 50, color = "steelblue", fill = "steelblue") +
  geom_vline(xintercept = base_model.phono.int.aic - psc_model.ld.aic, color = "darkred", size = 2) +
  ggtitle(" ", subtitle = subtitle_psc_subset.ld)  +
  labs(x = expression(Delta[AIC]), y = ' ') +
  coord_cartesian(xlim = (c(0, 560)), ylim = (c(0, 60))) +
  facet_grid(. ~ r)
(1+sum(stats_subset$deltaAIC.psc_ld_subset >= (base_model.phono.int.aic - psc_model.ld.aic)))/(length(iters)+1)

grid.arrange(p1.fsc_subset.aic, p2.fsc_subset.aic, ncol = 2)
grid.arrange(p3.fsc_subset.aic, p4.fsc_subset.aic, ncol = 2)

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


pca_te_model.ortho <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + morph + PC1_te + PC2_te + PC3_te, data = df.ortho
)
pca_te_model.ortho.aic = AIC(pca_te_model.ortho)
summary(pca_te_model.ortho)$r.sq
summary(pca_te_model.ortho)
base.ortho.aic - pca_te_model.ortho.aic


pca_ld_model.ortho <- lm(
  formula = aoa ~ n_chars.z + concr.z + val.z + freq.z + morph + PC1_ld + PC2_ld + PC3_ld, data = df.ortho
)
pca_ld_model.ortho.aic = AIC(pca_ld_model.ortho)
summary(pca_ld_model.ortho)$r.sq
summary(pca_ld_model.ortho)
base.ortho.aic - pca_ld_model.ortho.aic


##### PCA FOR NEIGHBORHOOD MEASURES - PHONOLOGY #####
neighborhoods_te.phono.pca <- prcomp(df.phono[c("PSC_te.z", "pnd.z", "snd.z")], center = TRUE, scale. = TRUE)
neighborhoods_ld.phono.pca <- prcomp(df.phono[c("PSC_ld.z", "pnd.z", "snd.z")], center = TRUE, scale. = TRUE)

plot(cumsum(neighborhoods_te.phono.pca$sdev^2 / sum(neighborhoods_te.phono.pca$sdev^2)), type="b", col='red')
lines(cumsum(neighborhoods_ld.phono.pca$sdev^2 / sum(neighborhoods_ld.phono.pca$sdev^2)), type="b", col='steelblue')

corrplot(cor(neighborhoods_te.phono.pca$x[,1:3], df.phono[c("PSC_te.z", "pnd.z", "snd.z")]))
pc_phono = df.phono[c("PSC_ld.z", "pnd.z", "snd.z")]
colnames(pc_phono) <- c(expression(PSC[ld]), "PND", "SND")
corrplot(cor(neighborhoods_ld.phono.pca$x[,1:3], pc_phono), 
         title = "Correlation matrix: PCs",
         mar = c(2,2,2,2),
         method = 'number',
         tl.col = 'black',
         tl.srt = 30)

df.phono['PC1_te_phono'] = neighborhoods_te.phono.pca$x[,1]
df.phono['PC2_te_phono'] = neighborhoods_te.phono.pca$x[,2]
df.phono['PC3_te_phono'] = neighborhoods_te.phono.pca$x[,3]

df.phono['PC1_ld_phono'] = neighborhoods_ld.phono.pca$x[,1]
df.phono['PC2_ld_phono'] = neighborhoods_ld.phono.pca$x[,2]
df.phono['PC3_ld_phono'] = neighborhoods_ld.phono.pca$x[,3]


pca_te_model.phono <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + morph + PC1_te_phono + PC2_te_phono + PC3_te_phono, data = df.phono
)
pca_te_model.phono.aic = AIC(pca_te_model.phono)
summary(pca_te_model.phono)$r.sq
summary(pca_te_model.phono)
base.phono.aic - pca_te_model.phono.aic


pca_ld_model.phono <- lm(
  formula = aoa ~ n_phons.z + concr.z + val.z + freq.z + morph + PC1_ld_phono + PC2_ld_phono + PC3_ld_phono, data = df.phono
)
pca_ld_model.phono.aic = AIC(pca_ld_model.phono)
summary(pca_ld_model.phono)$r.sq
summary(pca_ld_model.phono)
base.phono.aic - pca_ld_model.phono.aic
