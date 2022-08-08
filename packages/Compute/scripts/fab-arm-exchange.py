#name: Fab Arm Exchange
#language: python
#tags: antibody-eng
#meta.department: BTDS
#meta.status: Upstream
#meta.group: Bioengineering process
#input: double meaConcentration = 60 {category:FAE Reaction Parameters; units:mM; caption:2-MEA Concentration} [2-MEA Concentration]
#input: double phOfFaeReaction = 7.5 {category:FAE Reaction Parameters; caption:pH of FAE Reaction} [pH of FAE Reaction]
#input: double faeReactionTime = 3 {category:FAE Reaction Parameters; units:hr; caption:FAE Reaction Time} [FAE Reaction Time]
#input: double ionicStrenght = 100 {category:FAE Reaction Parameters; units:mM} [Ionic Strenght]
#input: double targetIgGConc = 10.5 {category:FAE Reaction Parameters; units:mg/ml; caption:Target IgG Conc} [Target IgG Conc]
#input: double massOverload = 1.05 {category:FAE Reaction Parameters; caption:Mass Overload (mAb1/mAb2)} [Mass Overload (mAb1/mAb2)]
#input: double P1MW = 148882 {category:FAE Reaction Parameters; units:mol/g; caption:Parental 1 MW} [Parental 1 MW]
#input: double reactionVolume = 0.5 {category:FAE Reaction Parameters; units:L} [Reaction Volume]
#input: double partialPressureOfOxygen = 0.21 {category:FAE Reaction Parameters; caption:Partial Pressure of Oxygen} [Partial Pressure of Oxygen]
#input: double doInitial = 0.9 {category:FAE Reaction Parameters; caption:DO (Initial)} [DO (Initial)]
#input: double doBuffer = 0.9 {category:FAE Reaction Parameters; caption:DO (Buffer)} [DO (Buffer)]
#input: double oxygenTransferRate = 0.443 {category:FAE Reaction Parameters; units:1/hr; caption:Oxygen Transfer rate} [Oxygen Transfer rate]
#input: double pressure = 1 {category:FAE Reaction Parameters; units:Bar} [Pressure]
#input: double temperature = 25 {category:FAE Reaction Parameters; units:C} [Temperature]
#input: bool Lid = False {category:FAE Reaction Parameters; caption:Lid (False=Closed, True=Opened)} [Lid (False=Closed, True=Opened)]
#input: double headSpace = 0.5 {category:FAE Reaction Parameters; units:L} [Head Space]
#input: double dfVolume = 12 {category:DF 1 Parameters; units:DV; caption:DF Volume} [DF Volume]
#input: double filterateRateDF = 10 {category:DF 1 Parameters; units:ml/min; caption:Filterate Rate (DF)} [Filterate Rate (DF)]
#input: double phOfUfdfBuffer = 7.5 {category:DF 1 Parameters; caption:pH of UFDF Buffer} [pH of UFDF Buffer]
#input: double phOfSecondUfdfBuffer = 7.5 {category:DF 2 Parameters (if any); caption:pH of 2nd UFDF Buffer} [pH of 2nd UFDF Buffer]
#input: double secondDfVolume = 0 {category:DF 2 Parameters (if any); units:DV; caption:2nd DF Volume} [2nd DF Volume]
#input: double dfConcentration = 30 {category: UF Parameters; caption:DF Concentration} [DF Concentration]
#input: double filterateRateUf = 10 {category: UF Parameters; units:ml/min; caption:Filterate Rate (UF)} [Filterate Rate (UF)]
#input: double holdTime = 0 {category:Hold Times Parameters (if any); units:hr} [Hold Time]
#input: double phDuringHoldTime = 7.5 {category:Hold Times Parameters (if any); caption:pH During Hold Time} [pH During Hold Time]
#input: double holdTime2 = 0 {category:Hold Times Parameters (if any); units:hr; caption:2nd Hold Time} [2nd Hold Time]
#input: double phDuringSecondHoldTime = 7.5 {category:Hold Times Parameters (if any); caption:pH During 2nd Hold Time} [pH During 2nd Hold Time]
#output: double c { category: Group A }
#output: double d { category: Group B }

c = holdTime + phOfSecondUfdfBuffer * massOverload
d = dfConcentration * phOfFaeReaction - meaConcentration