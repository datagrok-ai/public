/* eslint-disable max-len */
/** Bioreactor demo model */

import {LINK} from '../ui-constants';
import {Model} from '../model';

/** Bioreactor model specification */
const MODEL = `#name: Bioreactor
#description: Bioreactor simulation
#comment: 
  Source: https://doi.org/10.1074/jbc.RA117.000303.
#equations:

  d(FFox)/dt = -E11 + E12

  d(KKox)/dt = -E21 + E22

  d(FFred)/dt = E11 - E12 - E31 + E32

  d(KKred)/dt = E21 - E22 - E41 + E42

  d(Ffree)/dt = 2 * (E31 - E32) - E51 + E52

  d(Kfree)/dt = 2 * (E41 - E42) - E51 + E52

  d(FKred)/dt = E51 - E52 - E61 + E62

  d(FKox)/dt = E61 - E62

  d(MEAthiol)/dt = 2 * (-E11 + E12 - E21 + E22 + E31 + E41 - E32 - E42 - E62 - ktox * E71 * E72)
                   - (MEAthiol + MA) * (Fin + Fper) / VL

  d(CO2)/dt = (Fin * pO2sat * 0.0022 - 2 * Fper * CO2) / VL + OTR - 0.5 * ktox * E71 * E72

  d(yO2P)/dt = -OTR * (VL / (VTV - VL)) * R * T * P + yO2in * qin - yO2P * qout

  d(CYST)/dt = ktox * E71 * E72 - krcyst * CYST * E0 - (Fin + Fper) * CYST / VL

  d(VL)/dt = Fin - Fper

#expressions:

  KF = pow(VL, -0.65) * 0.065 * pow(speed**3 * diam**5 * power / 2.16e12, 0.361)

  MA = MEAthiol * pow(10, pH - pKa2MEA)

  qout = qin - KF * (yO2P * H - CO2) * VL * R * T / (P * 1000)
  
  OTR = KF * (yO2P * H - CO2)

  Fin = t < switchTime ? 0 : 0.025

  Fper = t < switchTime ? 0.025 : Fin

  E0 = VLinit / VL

  E1 = (MA * E0)**2
  
  E11 = k1red * FFox * E0 * E1

  E12 = k1ox * FFred * E0

  E31 = k2Fd * FFred * E0

  E32 = k2Fa * (Ffree * E0)**2 * E1

  E21 = k1red * KKox * E0 * E1

  E22 = k1ox * KKred * E0

  E41 = k2Kd * KKred * E0

  E42 = k2Ka * (Kfree * E0)**2 * E1

  E51 = k3FKa * Ffree * E0 * Kfree * E0

  E52 = k3FKd * FKred * E0

  E61 = k4ox * FKred * E0 * (CYST * E0)**2

  E62 = k4red * FKox * E0 * E1

  E70 = E0 * CO2

  E71 = (MEAthiol * E0)**2

  E72 = (E70 >= 0) ? sqrt(E70) : 0  

#argument: t, Reduction
  t0 = 0.0   {units: min; caption: Initial; category: Misc}                       [Initial time of simulation]
  t1 = 200   {units: min; caption: Reduction;   category: Duration; min: 200; max: 500}  [Time of reduction]
   h = 0.5   {units: min; caption: Step;    category: Misc; min: 0.1; max: 2}     [Time step of simulation]

#update:  Filtration
  duration = filtration

#inits:  
  FFox     = 0.2   {units: mmol/L; category: Initial values; min: 0.15; max: 0.25; step: 0.01}  [FF oxidized]
  KKox     = 0.2   {units: mmol/L; category: Initial values; min: 0.15; max: 0.25; step: 0.01}  [KK oxidized]
  FFred    = 0.1   {units: mmol/L; category: Initial values; min: 0.08; max: 0.12; step: 0.01}  [FF reduced]
  KKred    = 0.1   {units: mmol/L; category: Initial values; min: 0.08; max: 0.12; step: 0.01}  [KK reduced]
  Ffree    = 0     {units: mmol/L; category: Initial values}                                    [F free]
  Kfree    = 0     {units: mmol/L; category: Initial values}                                    [K free]
  FKred    = 0     {units: mmol/L; category: Initial values}                                    [FK reduced]
  FKox     = 0     {units: mmol/L; category: Initial values}                                    [FK oxidized]
  MEAthiol = 15    {units: mmol/L; category: Initial values; min: 10;   max: 16}                [MEAthiol]
  CO2      = 0.12  {units: mmol/L; category: Initial values; min: 0.09; max: 0.15}              [Dissolved oxygen]
  yO2P     = 0.209 {units: atm;    category: Initial values}                                    [Atm headspace]
  CYST     = 0     {units: mmol/L; category: Initial values}                                    [Cystamine]
  VL       = 7.2   {units: L;      category: Initial values}                                    [Liquid volume]

#output:
  t
  FFox
  KKox
  FFred
  KKred
  Ffree
  Kfree
  FKred
  FKox
  MEAthiol
  CO2
  yO2P
  CYST

#constants:
   VLinit = 7.2
      VTV = 10
    speed = 400
     diam = 6
    power = 2.1
       pH = 7.4
    k1red = 5.604e-2
     k1ox = 1.08e-2
     k2Fd = 1.35
     k2Fa = 1.104e8
     k2Kd = 4.038e-2
     k2Ka = 1.2e8
    k3FKa = 1.812e8
    k3FKd = 1.188e-2
     k4ox = 1.08e-2
    k4red = 5.604e-2
     ktox = 5e-3
   krcyst = 0
   pO2sat = 100
  pKa2MEA = 8.18
        H = 1.072069378
        R = 8.2e-2

#parameters:
  filtration = 300   {min: 100; max: 500; units: min; category: Duration} [Time of filtration]          
         qin =    1  {units: L/min; caption: Gas;         category: Parameters;  min: 0.5; max: 1.5}            [Gas to headspace]
       yO2in = 0.21  {              caption: O2 fraction; category: Parameters;  min: 0.1; max: 0.9}            [Oxygen mole fraction]
           T =  300  {units: K;     caption: temperature; category: Parameters;  min: 250; max: 350}            [System temperature]
           P =    1  {units: atm;   caption: pressure;    category: Parameters;  min: 1;   max: 2}              [Headspace pressure]
  switchTime =  135  {units: min;   caption: switch;   category: Misc;        min: 70;  max: 180; step: 10}  [Switch mode time]
  
#meta.inputs: mode {caption: Process mode; category: Process parameters; choices: OpenFile("System:AppData/DiffStudio/library/bioreactor-inputs.csv")} [Reactions flow mode]
`;

/** UI options */
const UI_OPTS = {
  inputsTabDockRatio: 0.17,
  graphsDockRatio: 0.85,
};

/** Info for help panel */
const INFO = `# Model
Simulation of a controlled fab-arm exchange kinetic
[mechanism](https://doi.org/10.1074/jbc.RA117.000303).
# Try
Interactive results based on input changes.
# Performance
1000 times faster than the previous version.
# Complexity
Each time you change inputs, a system of 12 nonlinear ordinary differential equations is solved.
# No-code
[Diff Studio](${LINK.DIF_STUDIO})
enables the creation of complex models without writing code.
# Learn more
* [Sensitivity analysis](${LINK.SENS_AN})
* [Parameter optimization](${LINK.FITTING})`;

/** Bioreactor demo model */
export const BIOREACTOR_DEMO = new Model(MODEL, UI_OPTS, INFO);
