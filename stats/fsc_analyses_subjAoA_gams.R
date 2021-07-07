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
library(mgcv)
library(GGally)
library(itsadug)
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

##### DATA PRE-PROCESSING #####
#load the data
data <- read.csv(file = "output/fsc_measures.csv")
iconicity <- read.csv(file = '/media/gioca90/University/TiU/Research/Resources/iconicity_ratings.csv')
df.mald = read.csv(file = '/media/gioca90/University/TiU/Research/Resources/MALD/AllData.txt', sep = '\t')

#transform morph into factor
data$morph <- factor(data$morph,    
                     level=c(0,1),    
                     labels=c("Mono","Poly"))


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


##### GAMs - ORTHOGRAPHY #####
# Base Model
base_model.ortho <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(old.z) + s(snd.z) + morph, data = df.ortho)
summary(base_model.ortho)$r.sq
# 0.5416503
summary(base_model.ortho)
# all smooths are significant
base.ortho.aic <- AIC(base_model.ortho)
# 34538.17
concurvity(base_model.ortho)
# some concurvity involving length and old


par(mfrow=c(2,3))
plot_smooth(base_model.ortho, view="freq.z", rug=TRUE, col="steelblue", main = "Orthography predictors (base model)",
            xlab = 'frequency', ylab = 'AoA')
plot_smooth(base_model.ortho, view="old.z", rug=TRUE, col="steelblue", 
            xlab = 'OLD20', ylab = 'AoA')
plot_smooth(base_model.ortho, view="n_chars.z", rug=TRUE, col="steelblue", 
            xlab = 'Length', ylab = 'AoA')
plot_smooth(base_model.ortho, view="concr.z", rug=TRUE, col="steelblue", 
            xlab = 'Concreteness', ylab = 'AoA')
plot_smooth(base_model.ortho, view="val.z", rug=TRUE, col="steelblue", 
            xlab = 'Valence', ylab = 'AoA')
plot_smooth(base_model.ortho, view="snd.z", rug=TRUE, col="steelblue", 
            xlab = 'SND', ylab = 'AoA')


# Base Model with snd*old interaction
base_model.ortho.int <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + ti(snd.z,old.z) + morph, data = df.ortho)
base_model.ortho.int.aic = AIC(base_model.ortho.int)
# 34536.03
summary(base_model.ortho.int)$r.sq
# 0.5419129
summary(base_model.ortho.int)
# tensor product not significant
base.ortho.aic - base_model.ortho.int.aic
# 2.140394

par(mfrow=c(1, 1))
vis.gam(base_model.ortho.int, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'old.z'), main = "composite smooth: OLD20, SND", xlab = "SND", ylab = "OLD20")


# Iconicity control
icon_model.ortho <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + morph + s(iconicity.z), data = df.ortho.icon)
icon_model.ortho.aic = AIC(icon_model.ortho)
# 7433.133
summary(icon_model.ortho)$r.sq
# 0.5912329
summary(icon_model.ortho)
# all smooths are significant except for length

plot_smooth(icon_model.ortho, view="iconicity.z", rug=TRUE, col="steelblue", main='Effect of iconicity on AoA (ortho)',
            xlab = 'Iconicity', ylab = 'AoA')


# target embedded neighbors
osc_model.te <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + morph + s(OSC_te.z), data = df.ortho)
osc_model.te.aic = AIC(osc_model.te)
# 34476.67
summary(osc_model.te)$r.sq
# 0.5452334
summary(osc_model.te)
# smooth of OSC_te is significant
base_model.ortho.int.aic - osc_model.te.aic
# 59.36626

plot_smooth(osc_model.te, view="OSC_te.z", rug=TRUE, col="steelblue", main='Effect of OSC[te] on AoA', xlab = expression('OSC'[te]), ylab = 'AoA')


# target embedded neighbors over iconicity
osc_icon_model.te <- bam(
  aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + morph + s(OSC_te.z) + s(iconicity.z), data = df.ortho.icon
  )
osc_icon_model.te.aic = AIC(osc_icon_model.te)
# 7429.523
summary(osc_icon_model.te)$r.sq
# 0.5928225
summary(osc_icon_model.te)
# both smooths for iconicity and OSCte are significant, although OSCte only slightly
icon_model.ortho.aic - osc_icon_model.te.aic
# 3.609972

par(mfrow=c(1, 2))
plot_smooth(osc_icon_model.te, view="OSC_te.z", rug=TRUE, col="steelblue", main='Effect of OSC[te] and Iconicity on AoA', 
            xlab = expression('OSC'[te]), ylab = 'AoA')
plot_smooth(osc_icon_model.te, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')


# OSCte*morph interaction
osc_by_morph_model.te <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + s(OSC_te.z, by = morph), data = df.ortho)
osc_by_morph_model.te.aic = AIC(osc_by_morph_model.te)
# 34568.79
summary(osc_by_morph_model.te)$r.sq
# 0.5405926
summary(osc_by_morph_model.te)
# OSCte smooth is significant for both morphological statuses, but AIC is worse
osc_model.te.aic - osc_by_morph_model.te.aic
# -92.12236

par(mfrow=c(1, 1))
plot_smooth(osc_by_morph_model.te, view="OSC_te.z", rug=TRUE, main='Interaction between OSC[te] and morphological status', 
            xlab = expression('OSC'[te]), ylab = 'AoA', plot_all='morph')



# Levenshtein neighbors
osc_model.ld <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + s(OSC_ld.z) + morph, data = df.ortho)
osc_model.ld.aic = AIC(osc_model.ld)
# 34252.9
summary(osc_model.ld)$r.sq
# 0.5569608
summary(osc_model.ld)
# smooth for OSCld is significant
base_model.ortho.int.aic - osc_model.ld.aic
# 283.1305

par(mfrow=c(1, 1))
plot_smooth(osc_model.ld, view="OSC_ld.z", rug=TRUE, col="steelblue", main='Effect of OSC[ld] on AoA', xlab = expression('OSC'[ld]), ylab = 'AoA')


# Levenshtein neighbors over iconicity
osc_icon_model.ld <- bam(
  aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + s(OSC_ld.z) + s(iconicity.z) + morph, data = df.ortho.icon
  )
osc_icon_model.ld.aic = AIC(osc_icon_model.ld)
# 7330.609
summary(osc_icon_model.ld)$r.sq
# 0.6130958
summary(osc_icon_model.ld)
# both smooths for OSCld and iconicity are significant
icon_model.ortho.aic - osc_icon_model.ld.aic
# 102.5239

par(mfrow=c(1, 2))
plot_smooth(osc_icon_model.ld, view="OSC_ld.z", rug=TRUE, col="steelblue", 
            main='Effect of OSC[ld] and Iconicity on AoA', xlab = expression('OSC'[ld]), ylab = 'AoA')
plot_smooth(osc_icon_model.ld, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')


# OSCld*morph interaction
osc_by_morph_model.ld <- bam(aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(snd.z) + s(old.z) + s(OSC_ld.z, by=morph), data = df.ortho)
osc_by_morph_model.ld.aic = AIC(osc_by_morph_model.ld)
# 34314.4
summary(osc_by_morph_model.ld)$r.sq
# 0.5536563
summary(osc_by_morph_model.ld)
# OSCte smooth is significant for both morphological statuses, but AIC is worse
osc_model.ld.aic - osc_by_morph_model.ld.aic
# -61.49383

par(mfrow=c(1, 1))
plot_smooth(osc_by_morph_model.ld, view="OSC_ld.z", rug=TRUE, xlab = expression('OSC'[ld]), ylab = 'AoA', plot_all='morph')



##### GAMs - PHONOLOGY #####
# Base Model
base_model.phono <- bam(aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph, data = df.phono)
summary(base_model.phono)$r.sq
# 0.5628369
summary(base_model.phono)
# all smooths are significant
base.phono.aic <- AIC(base_model.phono)
# 27236.46

par(mfrow=c(2,3))
plot_smooth(base_model.phono, view="freq.z", rug=TRUE, col="steelblue", main = "Phonology predictors (base model)",
            xlab = 'frequency', ylab = 'AoA')
plot_smooth(base_model.phono, view="pnd.z", rug=TRUE, col="steelblue", 
            xlab = 'PND', ylab = 'AoA')
plot_smooth(base_model.phono, view="n_phons.z", rug=TRUE, col="steelblue", 
            xlab = 'Length', ylab = 'AoA')
plot_smooth(base_model.phono, view="concr.z", rug=TRUE, col="steelblue", 
            xlab = 'Concreteness', ylab = 'AoA')
plot_smooth(base_model.phono, view="val.z", rug=TRUE, col="steelblue", 
            xlab = 'Valence', ylab = 'AoA')
plot_smooth(base_model.phono, view="snd.z", rug=TRUE, col="steelblue", 
            xlab = 'SND', ylab = 'AoA')


# Base Model with pnd*snd interaction
base_model.phono.int <- bam(aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + ti(pnd.z, snd.z) + morph, data = df.phono)
summary(base_model.phono.int)$r.sq
# 0.5631592
summary(base_model.phono.int)
# tensor product not significant
base_model.phono.int.aic <- AIC(base_model.phono.int)
# 27237.88
base.phono.aic - base_model.phono.int.aic
# -1.427511

par(mfrow=c(1, 1))
vis.gam(base_model.phono.int, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'pnd.z'), main = "composite smooth: PND, SND", xlab = "SND", ylab = "PND")


# Iconicity control
icon_model.phono <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z), data = df.phono.icon
  )
icon_model.phono.aic = AIC(icon_model.phono)
# 6516.004
summary(icon_model.phono)$r.sq
# 0.591665
summary(icon_model.phono)
# iconicity smooth is significant

plot_smooth(icon_model.phono, view="iconicity.z", rug=TRUE, col="steelblue", main='Effect of iconicity on AoA (phono)',
            xlab = 'Iconicity', ylab = 'AoA')


# target embedded neighbors
psc_model.te <- bam(aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(PSC_te.z), data = df.phono)
psc_model.te.aic = AIC(psc_model.te)
# 27184.72
summary(psc_model.te)$r.sq
# 0.566543
summary(psc_model.te)
# smooth for PSC_te is significant
base_model.phono.int.aic - psc_model.te.aic
# 53.15989

plot_smooth(psc_model.te, view="PSC_te.z", rug=TRUE, col="steelblue", main='Effect of PSC[te] on AoA', xlab = expression('PSC'[te]), ylab = 'AoA')


# target embedded neighbors over iconicity
psc_icon_model.te <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z) + s(PSC_te.z), 
  data = df.phono.icon
  )
psc_icon_model.te.aic = AIC(psc_icon_model.te)
# 6518.418
summary(psc_icon_model.te)$r.sq
# 0.5918991
summary(psc_icon_model.te)
# smooth for PSCte no longer significant after including iconicity
icon_model.phono.aic - psc_icon_model.te.aic
# -2.41442

par(mfrow=c(1, 2))
plot_smooth(psc_icon_model.te, view="PSC_te.z", rug=TRUE, col="steelblue", 
            main='Effect of PSC[te] abd Iconicity on AoA', xlab = expression('PSC'[te]), ylab = 'AoA')
plot_smooth(psc_icon_model.te, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')


# PSCte*morph interaction
psc_by_morph_model.te <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + s(PSC_te.z, by=morph), data = df.phono
  )
psc_by_morph_model.te.aic = AIC(psc_by_morph_model.te)
# 27297.6
summary(psc_by_morph_model.te)$r.sq
# 0.5591715
summary(psc_by_morph_model.te)
# PSCte smooths are significant regardless of morphological status, but worse AIC
psc_model.te.aic - psc_by_morph_model.te.aic
# -112.8751

par(mfrow=c(1, 1))
plot_smooth(psc_by_morph_model.te, view="PSC_te.z", rug=TRUE, xlab = expression('PSC'[te]), ylab = 'AoA', plot_all='morph')



# Levenshtein neighbors
psc_model.ld <- bam(aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(PSC_ld.z), data = df.phono)
psc_model.ld.aic = AIC(psc_model.ld)
# 26817.51
summary(psc_model.ld)$r.sq
# 0.5892933
summary(psc_model.ld)
# smooth for PSCld is significant
base_model.phono.int.aic - psc_model.ld.aic
# 420.3715

plot_smooth(psc_model.ld, view="PSC_ld.z", rug=TRUE, col="steelblue", main='Effect of PSC[ld] on AoA', xlab = expression('PSC'[ld]), ylab = 'AoA')


# Levenshtein neighbors over iconicity
psc_icon_model.ld <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(iconicity.z) + s(PSC_ld.z), 
  data = df.phono.icon
  )
psc_icon_model.ld.aic = AIC(psc_icon_model.ld)
# 6422.73
summary(psc_icon_model.ld)$r.sq
# 0.6137275
summary(psc_icon_model.ld)
# both smooths are significant
icon_model.phono.aic - psc_icon_model.ld.aic
# 93.27363

par(mfrow=c(1, 2))
plot_smooth(psc_icon_model.ld, view="PSC_ld.z", rug=TRUE, col="steelblue", 
            main='Effect of PSC[ld] abd Iconicity on AoA', xlab = expression('PSC'[ld]), ylab = 'AoA')
plot_smooth(psc_icon_model.ld, view="iconicity.z", rug=TRUE, col="steelblue", xlab = 'Iconicity', ylab = 'AoA')


# PSCld*morph interaction
psc_by_morph_model.ld <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + s(PSC_ld.z, by=morph), data = df.phono
  )
psc_by_morph_model.ld.aic = AIC(psc_by_morph_model.ld)
# 26815.48
summary(psc_by_morph_model.ld)$r.sq
# 0.5899462
summary(psc_by_morph_model.ld)
# PSCld smooths are significant regardless of morphological status, but only slightly better AIC due to overfitting at high PSCld for polymorphemic words
psc_model.ld.aic - psc_by_morph_model.ld.aic
# 2.034638

par(mfrow=c(1, 1))
plot_smooth(psc_by_morph_model.ld, view="PSC_ld.z", rug=TRUE, xlab = 'PSC_ld', ylab = 'AoA', plot_all='morph', ylim = c(6, 15))



##### Interaction FSC*SND #####

# 2-way interaction: SND*OSC
osc_by_snd_model.ld <- bam(
  aoa ~ s(n_chars.z) + s(concr.z) + s(val.z) + s(freq.z) + s(old.z) + s(snd.z) + morph + s(OSC_ld.z) + ti(OSC_ld.z,snd.z), data = df.ortho
  )
osc_by_snd_model.ld.aic = AIC(osc_by_snd_model.ld)
# 34255.5
summary(osc_by_snd_model.ld)$r.sq
# 0.5570917
summary(osc_by_snd_model.ld)
# tensor product between OSCld and snd not significant
osc_model.ld.aic - osc_by_snd_model.ld.aic
# -2.600942

par(mfrow=c(1, 1))
vis.gam(osc_by_snd_model.ld, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'OSC_ld.z'), main = "composite smooth: SND, OSC_ld", xlab = "SND", ylab = "OSC_ld")


# 2-way interaction: SND*PSC
psc_by_snd_model.ld <- bam(
  aoa ~ s(n_phons.z) + s(concr.z) + s(val.z) + s(freq.z) + s(pnd.z) + s(snd.z) + morph + s(PSC_ld.z) + ti(PSC_ld.z,snd.z), data = df.phono
  )
psc_by_snd_model.ld.aic = AIC(psc_by_snd_model.ld)
# 26817.94
summary(psc_by_snd_model.ld)$r.sq
# 0.589291
summary(psc_by_snd_model.ld)
# tensor product between PSCld and snd not significant
psc_model.ld.aic - psc_by_snd_model.ld.aic
# -0.4306508

par(mfrow=c(1, 1))
vis.gam(psc_by_snd_model.ld, too.far=0.1, color="topo", plot.type="contour", 
        view=c('snd.z', 'PSC_ld.z'), main = "composite smooth: SND, PSC_ld", xlab = "SND", ylab = "PSC_ld")
