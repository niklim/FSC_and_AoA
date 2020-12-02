#load packages
install.packages("ranger")
install.packages("effects")
install.packages("MASS")
install.packages("psych")
install.packages("lsr")
install.packages("dplyr")
install.packages("car")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("rockchalk")
install.packages("ggeffects")
install.packages("plyr")
install.packages("vip")

library(plyr)
library(vip)
library(ranger)
library(effects)
library(ggeffects)
library(MASS)
library(psych)
library(lsr)
library(dplyr)
library(car)
library(reshape2)
library(ggplot2)
library(rockchalk)
set.seed(1)

#function to transform all variables into normal shape 
box_cox_transf <- function(x, min, max) {
  bc = boxcox(x ~ 1, lambda = seq(min, max, 0.1))
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

#load the data
data <- read.csv(file = "/Users/niklaslimacher/Documents/Thesis/Code/data_morph.csv")

#remove first column
data <- select(data, -X)

#transform morph into factor
data$morph <- factor(data$morph,    
                     level=c(0,1),    
                     labels=c("Mono","Poly"))

#get overview of the variables
head(data, 5)
tail(data, 5)
summary(data)
describe(data)


### Transform variables into normal shape and comparable units (SD) ###
df <- data[,FALSE]

df$aoa = data$aoa
df$OSC_m.z = c(scale(box_cox_transf(data$OSC_m + 0.0001, -10, 10), center = T, scale = T))
df$OSC_h.z = c(scale(box_cox_transf(data$OSC_h, -10, 10), center = T, scale = T))
df$PSC_m.z = c(scale(box_cox_transf(data$PSC_m + 0.0001, -10, 10), center = T, scale = T))
df$PSC_h.z = c(scale(box_cox_transf(data$PSC_h, -10, 10), center = T, scale = T))
df$conc.z = c(scale(box_cox_transf(data$conc, -10, 10), center = T, scale = T))
df$val.z = c(scale(box_cox_transf(data$val, -10, 10), center = T, scale = T))
df$wordlength.z = c(scale(box_cox_transf(data$word.length, -10, 10), center = T, scale = T))
df$freq.z = c(scale(box_cox_transf(data$freq_subtlex, -10, 10), center = T, scale = T))
df$old.z = c(scale(box_cox_transf(data$old20, -10, 10), center = T, scale = T))
df$morph = data$morph

#get histograms of variabbles
hist(data$aoa)

hist(data$word.length)
hist(df$wordlength.z)

hist(data$OSC_m, breaks = 30, main = "Histogram of OSC (MA) Values", xlab = "OSC Value")
hist(df$OSC_m.z, breaks = 30, main = "Histogram of Transformed OSC (MA) Values", xlab = "OSC Values")

?hist

hist(data$OSC_h, breaks = 30)
hist(df$OSC_h.z, breaks = 30)

hist(data$PSC_m, breaks = 30, main = "Histogram of PSC (MA) Values", xlab = "OSC Value")
hist(df$PSC_m.z, breaks = 30, main = "Histogram of Transformed PSC (MA) Values", xlab = "OSC Value")

hist(data$PSC_h, breaks = 30)
hist(df$PSC_h.z, breaks = 30)

hist(data$conc)
hist(df$conc.z)

hist(data$val)
hist(df$val.z)

hist(data$freq_subtlex)
hist(df$freq.z)

mono <- filter(df, morph == "Mono")
poly <- filter(df, morph == "Poly")

hist(mono$OSC_m.z, breaks = 30, main = "Histogram of OSC (MA) Values for Monomorphemic Words", xlab = "OSC Values" )
hist(poly$OSC_m.z, breaks = 30, main = "Histogram of OSC (MA) Values for Polymorphemic Words", xlab = "OSC Values")

hist(mono$PSC_m.z, breaks = 30, main = "Histogram of PSC (MA) Values for Monomorphemic Words", xlab = "PSC Values" )
hist(poly$PSC_m.z, breaks = 30, main = "Histogram of PSC (MA) Values for Polymorphemic Words", xlab = "PSC Values")

#check correlations & compute matrix
correlate(data, test = TRUE)
cor_mat <- round(cor(df[,-11]),3)
cor_mat
upper_tri <- get_upper_tri(cor_mat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
levels(melted_cormat$Var1)
melted_cormat$Var1 = plyr::mapvalues(melted_cormat$Var1, 
                                     from = c("aoa","OSC_m.z", "OSC_h.z","PSC_m.z","PSC_h.z","conc.z","val.z","wordlength.z","freq.z","old.z"), 
                                     to = c("AoA", "OSC_m", "OSC_h", "PSC_m", "PSC_h", "Concr","Val","Len","Freq", "OLD20" ))
melted_cormat$Var2 = plyr::mapvalues(melted_cormat$Var2, 
                                     from = c("aoa","OSC_m.z", "OSC_h.z","PSC_m.z","PSC_h.z","conc.z","val.z","wordlength.z","freq.z","old.z"), 
                                     to = c("AoA", "OSC_m", "OSC_h", "PSC_m", "PSC_h", "Concr","Val","Len","Freq", "OLD20" ))

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(high = "darkred", low = "steelblue", mid = "lightgray", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 11, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 4) +
  ggtitle("Pairwise Correlations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.75),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 10, barheight = 2,
                               title.position = "top", title.hjust = 0.5))



a <- ggplot(data, aes(aoa,OSC_h))
a+geom_point(shape = 20, color = "#9A3C3C") +
  ggtitle("AoA vs. OSC (H)")+
  xlab("AoA")+
  ylab("OSC (H)")

b <- ggplot(data, aes(aoa,OSC_m))
b + geom_point(shape = 20, color = "#9A3C3C") +
  ggtitle("AoA vs. OSC (M)")+
  xlab("AoA")+
  ylab("OSC (M)")

#################### Multiple Linear Regression Models ################################


### Base Model ###
base_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph, data = df)
AIC(base_model)
summary(base_model)$r.sq
summary(base_model)
vif(base_model)
base.AIC <- AIC(base_model)

### OSC Marelli ###

osc_m_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_m.z, data = df)
AIC(osc_m_model)
summary(osc_m_model)$r.sq
summary(osc_m_model)
vif(osc_m_model)
base.AIC - AIC(osc_m_model)


### OSC Hendrix ###
osc_h_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_h.z, data = df)
AIC(osc_h_model)
summary(osc_h_model)$r.sq
summary(osc_h_model)
base.AIC - AIC(osc_h_model)

### PSC Marelli ### 
psc_m_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + PSC_m.z, data = df)
AIC(psc_m_model)
summary(psc_m_model)$r.sq
summary(psc_m_model)
base.AIC - AIC(psc_m_model)


### PSC Hendrix ###
psc_h_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + PSC_h.z, data = df)
AIC(psc_h_model)
summary(psc_h_model)$r.sq
summary(psc_h_model)
base.AIC - AIC(psc_h_model)


### OSC & PSC Marelli ### 
opsc_m_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_m.z + PSC_m.z, data = df)
AIC(opsc_m_model)
summary(opsc_m_model)$r.sq
summary(opsc_m_model)
AIC(osc_m_model) - AIC(opsc_m_model)
AIC(psc_m_model) - AIC(opsc_m_model)

### OSC & PSC Hendrix ### 
opsc_h_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_h.z + PSC_h.z, data = df)
AIC(opsc_h_model)
summary(opsc_h_model)$r.sq
summary(opsc_h_model)
AIC(osc_h_model) - AIC(opsc_h_model)
AIC(psc_h_model) - AIC(opsc_h_model)

#################### Interaction Plots OSC x Morph ################################

### Interaction OSC Marelli x Morph ###
oscxmorph_m_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + OSC_m.z * morph, data = df)
AIC(oscxmorph_m_model)
summary(oscxmorph_m_model)$r.sq
summary(oscxmorph_m_model)
AIC(osc_m_model) - AIC(oscxmorph_m_model)

lm.oscxmorph_m = data.frame(ggpredict(oscxmorph_m_model, terms = c("OSC_m.z","morph")))

ggplot(
  lm.oscxmorph_m, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line( 
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('OSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ OSC"[MA]*" x Morphological Complexity")) +
  labs(colour = "Morphological Complexity", fill = "Morphological Complexity")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

### Interaction OSC Hendrix x Morphological Complexity ###
oscxmorph_h_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + OSC_h.z * morph, data = df)
AIC(oscxmorph_h_model)
summary(oscxmorph_h_model)$r.sq
summary(oscxmorph_h_model)
AIC(osc_h_model) - AIC(oscxmorph_h_model)

lm.oscxmorph_h = data.frame(ggpredict(oscxmorph_h_model, terms = c("OSC_h.z","morph")))

ggplot(
  lm.oscxmorph_h, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line( 
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('OSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ OSC"[HS]*" x Morphological Complexity")) +
  labs(colour = "Morphological Complexity", fill = "Morphological Complexity")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


### Interaction PSC Marelli x Morph ###
pscxmorph_m_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + PSC_m.z * morph, data = df)
AIC(pscxmorph_m_model)
summary(pscxmorph_m_model)$r.sq
summary(pscxmorph_m_model)
AIC(psc_m_model) - AIC(pscxmorph_m_model)

lm.pscxmorph_m = data.frame(ggpredict(pscxmorph_m_model, terms = c("PSC_m.z","morph")))

ggplot(
  lm.pscxmorph_m, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line( 
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('PSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ PSC"[MA]*" x Morphological Complexity")) +
  labs(colour = "Morphological Complexity", fill = "Morphological Complexity")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

### Interaction OSC Hendrix x Morphological Complexity ###
pscxmorph_h_model <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + PSC_h.z * morph, data = df)
AIC(pscxmorph_h_model)
summary(pscxmorph_h_model)$r.sq
summary(pscxmorph_h_model)
AIC(psc_h_model) - AIC(pscxmorph_h_model)

lm.pscxmorph_h = data.frame(ggpredict(pscxmorph_h_model, terms = c("PSC_h.z","morph")))

ggplot(
  lm.pscxmorph_h, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line( 
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('PSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ PSC"[HS]*" x Morphological Complexity")) +
  labs(colour = "Morphological Complexity", fill = "Morphological Complexity")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 




#################### Interaction Plots Frequency x O/PSC ################################

### Model Interaction Frequency & OSC Marelli ###
m1 <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + old.z + morph + freq.z * OSC_m.z, data = df)
summary(m1)
AIC(m1)

lm.freqxosc_m = data.frame(ggpredict(m1, terms = c("OSC_m.z","freq.z[-1.476,-0.082,1.721")))
lm.freqxosc_m$group <- mapvalues(lm.freqxosc_m$group, from = c("-1", "0", "1"), to = c("Low", "Mid", "High"))
quantile(df$freq.z, c(.05,0.5,0.95))

ggplot(
  lm.freqxosc_m, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line(
     alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('OSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ Frequency x OSC"[MA])) +
  labs(colour = "Frequency", fill = "Frequency")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

### Interaction Frequency x OSC Hendrix ###
m2 <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + old.z + morph + freq.z * OSC_h.z, data = df)
summary(m2)
AIC(m2)

lm.freqxosc_h = data.frame(ggpredict(m2, terms = c("OSC_h.z","freq.z[-1.476,-0.082,1.721")))
lm.freqxosc_h$group <- mapvalues(lm.freqxosc_h$group, from = c("-1", "0", "1"), to = c("Low", "Mid", "High"))

ggplot(
  lm.freqxosc_h, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line( 
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('OSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ Frequency x OSC"[HS])) +
  labs(colour = "Frequency", fill = "Frequency")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


### Model Interaction Frequency & PSC Marelli ###
m3 <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + old.z + morph + freq.z * PSC_m.z, data = df)
summary(m3)
AIC(m3)

lm.freqxpsc_m = data.frame(ggpredict(m3, terms = c("PSC_m.z","freq.z[-1.476,-0.082,1.721")))
lm.freqxpsc_m$group <- mapvalues(lm.freqxpsc_m$group, from = c("-1", "0", "1"), to = c("Low", "Mid", "High"))

ggplot(
  lm.freqxpsc_m, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line(
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('PSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ Frequency x PSC"[MA])) +
  labs(colour = "Frequency", fill = "Frequency")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


### Model Interaction Frequency & PSC Hendrix ###
m4 <- lm(formula = aoa ~ wordlength.z + conc.z + val.z + old.z + morph + freq.z * PSC_h.z, data = df)
summary(m4)
AIC(m4)

lm.freqxpsc_h = data.frame(ggpredict(m4, terms = c("PSC_h.z","freq.z[-1.476,-0.082,1.721")))
lm.freqxpsc_h$group <- mapvalues(lm.freqxpsc_h$group, from = c("-1", "0", "1"), to = c("Low", "Mid", "High"))

ggplot(
  lm.freqxpsc_h, 
  aes(x = x, y = predicted, colour = group)
) +  
  geom_line(
    alpha = 0.9
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, 
        fill = group), alpha = 0.25
  ) +
  xlab('PSC') +
  ylab('AoA') +
  ggtitle(expression("AoA ~ Frequency x PSC"[HS])) +
  labs(colour = "Frequency", fill = "Frequency")+
  theme(
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 


#################### Interaction Effects PSC x Morph ################################

plotSlopes(pscxmorph_m_model, plotx = "PSC_m.z", modx = "morph", plotPoints = FALSE)

plotSlopes(pscxmorph_h_model, plotx = "PSC_h.z", modx = "morph", plotPoints = FALSE)

#################### Random Forest Regression Models ################################

### Random Forest Base ###
rf_base <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + morph + old.z, data = df, mtry = 5/3, importance = "impurity" )
print(rf_base)
vip(rf_base, aesthetics = list(col = "darkblue", fill = "darkblue"))

?ranger
### Random Forest OSC Marelli ###
rf_OSC_m <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_m.z, data = df, mtry = 6/3, importance = "impurity" )
print(rf_OSC_m)
vip(rf_OSC_m, aesthetics = list(col = "darkblue", fill = "darkblue"))


### Random Forest OSC Hendrix ###
rf_OSC_h <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_h.z, data = df, mtry = 6/3, importance = "impurity" )
print(rf_OSC_h)
vip(rf_OSC_h, aesthetics = list(col = "darkblue", fill = "darkblue"))


### Random Forest PSC Marelli ###
rf_PSC_m <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + PSC_m.z, data = df, mtry = 6/3, importance = "impurity" )
print(rf_PSC_m)
vip(rf_PSC_m, aesthetics = list(col = "darkblue", fill = "darkblue"))

### Random Forest PSC Hendrix ###
rf_PSC_h <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph +  PSC_h.z, data = df, mtry = 6/3, importance = "impurity" )
print(rf_PSC_h)
vip(rf_PSC_h, aesthetics = list(col = "darkblue", fill = "darkblue"))


### Random Forest OSC & PSC Marelli ###
rf_OPSC_m <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_m.z + PSC_m.z, data = df, mtry = 7/3, importance = "impurity" )
print(rf_OPSC_m)
vip(rf_OPSC_m, aesthetics = list(col = "lightblue3", fill = "lightblue3")) +
  ggtitle("Random Forest Variable Importance (MA)") +
  theme(axis.text=element_text(size=12))


### Random Forest OSC & PSC Hendrix ###
rf_OPSC_h <- ranger(formula = aoa ~ wordlength.z + conc.z + val.z + freq.z + old.z + morph + OSC_h.z + PSC_h.z, data = df, mtry = 7/3, importance = "impurity" )
print(rf_OPSC_h)
vip(rf_OPSC_h, aesthetics = list(col = "lightblue3", fill = "lightblue3")) +
  ggtitle("Random Forest Variable Importance (HS)") + 
  theme(axis.text=element_text(size=12))

citation()
