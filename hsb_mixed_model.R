library(lme4)
library(mlmhelpr)

data("hsb", package = "mlmhelpr")
head(hsb)

baseline_model <- lmer(mathach ~ ses + minority + female + catholic + size + 
                         pracad + disclim + himinty + meanses +
                         (1 + ses + minority + female | id), 
                       data = hsb)
summary(baseline_model)
