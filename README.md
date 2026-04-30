# Masterarbeit: Variablenselektion in Linearen Gemischten Modellen

Verglichen werden drei Ansätze:
- Bondell's Joint Selection Method
- Lin's Two-Stage Method
- Benchmark Method (basierend auf `glmmTMB`)

Die Implementierungen der Verfahren befinden sich im `R/`-Ordner:
- **Bondell's Methode** (Adaptive Group Lasso mit Cholesky-Zerlegung): `R/method_bondell.R`
- **Lin's Two-Stage Methode** (Penalized REML und penalized profile log-likelihood): `R/method_lin.R`
- **Jacob's Methode** (Benchmark-Modell mit L1-Penalty auf Fixed Effects): `R/fitModel.R` und `R/fitGlmmLasso.R`

Alle Methoden werden über den Wrapper `R/method_wrappers.R` aufgerufen.

Das Skript zur Ausführung der Simulation ist:
`scripts/main_simulation.R`

Ablauf der Simulation:
1. Generierung von Datensätzen (mithilfe von `R/generate_data.R`)
2. Anwendung der drei Methoden auf die generierten Daten
3. Berechnung der Performance-Metriken
4. Speicherung der Daten als CSV im `data/`
