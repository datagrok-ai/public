# Model
Simulation of a two-compartment pharmacokinetic-pharmacodynamic
([PK-PD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046)).
# Try
Interactive results based on input changes.
# Performance
Nonlinear systems of differential equations are solved within milliseconds.
# Quick reference

| Input | Meaning |
|---|---|
| Interval | Hours between doses |
| Dose | Amount per dose |
| Count | Number of doses |
| Central (concentration) | Central concentration at start |
| Peripheral (concentration) | Peripheral concentration at start |
| Rate constant | Drug absorption rate |
| Clearance | Drug elimination from blood |
| Central volume | Central compartment volume |
| Inter rate | Central-peripheral exchange rate |
| Peri volume | Peripheral compartment volume |
| Effect | Half-maximal effect concentration |
| Rate | Effect turnover rate |
| Begin | Dosing interval start |
| Step | ODE solver step size |
| Depo | Initial drug at absorption site |
| Central (amount) | Initial central compartment amount |
| Peripheral (amount) | Initial peripheral compartment amount |
| Init effect | Starting biomarker level |
