## Modell 1

# Daten generieren
{
  beta_1 <- c(5, 2)
  covX_1 <- matrix(1, nrow = 1) # p=1 -> x1
  covU_1 <- matrix(4)
  
  set.seed(123)
  data_1 <- generateData(
    beta = beta_1, n = 500, covX = covX_1, n_groups = 25, 
    grV = 1, slope = 0, covU = covU_1, sigma_eps = 1.5
  )
  
  print(head(data_1$data))
}

# Modell erstellen
{
  model_1 <- lmer(y ~ x1 + (1 | G1), data = data_1$data)
  print(summary(model_1))
}


## Modell 2
# Datensatz generieren
{
  beta_2 <- c(10, 3, -1)
  covX_2 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2) 
  covU_2 <- matrix(c(9, 1, 1, 4), nrow = 2) 
  
  set.seed(123)
  data_2 <- generateData(
    beta = beta_2, n = 1000, covX = covX_2, n_groups = 40, 
    grV = 2, slope = 1, covU = covU_2, sigma_eps = 2.0
  )
  
  print(head(data_2$data))
}

# Modell erstellen
{
  model_2 <- lmer(y ~ x1 + x2 + (1 + x1 | G1) + (1 | G2), data = data_2$data)
  print(summary(model_2))
}
