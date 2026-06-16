import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';
import {runDiffStudioModel} from './ivp-runtime';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: dock
export function dock() : void {
  PackageFunctions.dock();
}

//input: object problem 
//output: dataframe result
export function solve(problem: any) : any {
  return PackageFunctions.solve(problem);
}

//input: object problem 
//input: object options 
//output: dataframe result
export function solveEquations(problem: any, options: any) : any {
  return PackageFunctions.solveEquations(problem, options);
}

//name: DiffStudio Facet
//description: Faceted grid of line charts, one per output variable, for Diff Studio solutions
//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function diffStudioFacetViewer() : any {
  return PackageFunctions.diffStudioFacetViewer();
}

//name: Diff Studio
//description: Solver of ordinary differential equations systems
//output: view result
//meta.role: app
//meta.browsePath: Compute
export async function runDiffStudio() : Promise<any> {
  return await PackageFunctions.runDiffStudio();
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: runDiffStudioDemo() //wait: 100 
export async function runDiffStudioDemo() : Promise<void> {
  await PackageFunctions.runDiffStudioDemo();
}

//input: string content 
//meta.role: fileHandler
//meta.ext: ipv
export async function ivpFileHandler(content: string) : Promise<void> {
  await PackageFunctions.ivpFileHandler(content);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ivp
export async function previewIvp(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewIvp(file);
}

//input: dynamic treeNode 
export async function runDiffStudioTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.runDiffStudioTreeBrowser(treeNode);
}

//name: Ball flight
//description: Ball flight simulation
//tags: model
//input: double dB = 0.01 { category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3; minFormula: roB / 20000; maxFormula: roB / 4000 }
//input: double roB = 200 { category: Ball; caption: Density; description: Material density; units: kg/m^3; min: 200; max: 1200 }
//input: double v = 50 { category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec }
//input: double a = 45 { category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg }
//output: double maxDist { caption: Max distance }
//output: double maxHeight { caption: Max height }
//output: dataframe df { caption: Trajectory; viewer: Line chart(multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid() }
//editor: Compute2:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ball.png
//meta.help: ball-flight.md
//meta.dockSpawnConfig: {"Trajectory / Grid": {"dock-spawn-dock-ratio": 0.3, "dock-spawn-dock-type": "right", "dock-spawn-dock-to": "Trajectory / Line chart"}, "Output": {"dock-spawn-dock-ratio": 0.15, "dock-spawn-dock-type": "down", "dock-spawn-dock-to": "Trajectory / Line chart"}}
export function ballFlight(dB: number, roB: number, v: number, a: number) {
  return PackageFunctions.ballFlight(dB, roB, v, a);
}

//description: Return serialized initial value problem for ordinary differential equations
//input: string problem 
//output: object serialization
export function serializeEquations(problem: string) : any {
  return PackageFunctions.serializeEquations(problem);
}

//description: Perform ODEs serialization to JS-code
//input: dynamic serialization 
//output: string result
export function odesToCode(serialization: any) : string {
  return PackageFunctions.odesToCode(serialization);
}

//description: Solve initial value problem for ordinary differential equations
//input: string problem 
//output: dataframe result
export async function solveODE(problem: string) : Promise<any> {
  return await PackageFunctions.solveODE(problem);
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//output: dynamic result
//meta.demoPath: Compute | PK-PD Modeling
//test: demoSimPKPD() //wait: 100 
export async function demoSimPKPD() : Promise<any> {
  return await PackageFunctions.demoSimPKPD();
}

//name: Bioreactor Demo
//description: In-browser simulation of controlled fab-arm exchange mechanism
//output: dynamic result
//meta.demoPath: Compute | Bioreactor
//test: demoBioreactor() //wait: 100 
export async function demoBioreactor() : Promise<any> {
  return await PackageFunctions.demoBioreactor();
}

//description: Run model with Diff Studio UI
//input: string model 
//input: int inputsTabDockRatio 
//input: int graphsDockRatio 
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number) : Promise<void> {
  await PackageFunctions.runModel(model, inputsTabDockRatio, graphsDockRatio);
}

//name: Acid Production
//description: Gluconic acid (GA) production by Aspergillus niger modeling
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ga-production.png
//meta.help: library/acid-production.md
//meta.role: model
//meta.diffStudioModel: System:AppData/DiffStudio/library/acid-production.ivp
//editor: Compute2:RichFunctionViewEditor
//input: double _t0 = 0 {units: h; caption: initial; category: Misc} [Start of the process]
//input: double _t1 = 60 {units: h; caption: 1-st stage; category: Durations; min: 20; max: 80} [Duration of the 1-st stage]
//input: double _h = 0.1 {units: h; caption: step; category: Misc; min: 0.01; max: 1} [Time step of simulation]
//input: double X = 5 {units: kg/m³; caption: biomass; category: Initial concentrations; min: 1; max: 10} [Aspergillus niger biomass]
//input: double S = 150 {units: kg/m³; caption: glucose; category: Initial concentrations; min: 50; max: 200} [Glucose]
//input: double O = 7 {units: kg/m³; caption: oxygen; category: Initial concentrations; min: 1; max: 10} [Dissolved oxygen]
//input: double P = 0 {units: kg/m³; caption: acid; category: Initial concentrations; min: 0; max: 0.1} [Gluconic acid]
//input: double overall = 100 {units: h; category: Durations; min: 100; max: 140} [Overall duration]
//input: double muM = 0.668 {units: 1/h; category: Parameters} [Monod type model parameter]
//input: double alpha = 2.92 {category: Parameters} [Monod type model parameter]
//input: double beta = 0.131 {units: 1/h; category: Parameters} [Monod type model parameter]
//input: double gamma = 2.12 {category: Parameters} [Monod type model parameter]
//input: double lambda = 0.232 {units: 1/h; category: Parameters} [Monod type model parameter]
//input: double delta = 0.278 {category: Parameters} [Monod type model parameter]
//input: double phi = 0.00487 {units: 1/h; category: Parameters} [Monod type model parameter]
//input: double Ks = 130.9 {units: g/L; category: Parameters} [Monod type model parameter]
//input: double Ko = 0.000363 {units: g/L; category: Parameters} [Monod type model parameter]
//input: double Kla = 0.017 {units: 1/s; category: Parameters} [Volumetric mass transfer coefficient]
//input: double Cod = 15 {units: kg/m³; category: Parameters} [Liquid phase dissolved oxygen saturation concentration]
//output: dataframe df {caption: Acid Production; viewer: Grid(block: 100) | Line chart(block: 100, multiAxis: "true", segmentColumnName: "_Stage", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | DiffStudio Facet(block: 100, segmentColumnName: "_Stage")}
export async function ivpModel_Acid_Production(_t0: number, _t1: number, _h: number, X: number, S: number, O: number, P: number, overall: number, muM: number, alpha: number, beta: number, gamma: number, lambda: number, delta: number, phi: number, Ks: number, Ko: number, Kla: number, Cod: number): Promise<DG.DataFrame> {
  return await runDiffStudioModel('System:AppData/DiffStudio/library/acid-production.ivp', {_t0, _t1, _h, X, S, O, P, overall, muM, alpha, beta, gamma, lambda, delta, phi, Ks, Ko, Kla, Cod});
}

//name: Pollution
//description: The chemical reaction part of the air pollution model developed at The Dutch National Institute of Public Health and Environmental Protection
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/pollution.png
//meta.help: library/pollution.md
//meta.role: model
//meta.diffStudioModel: System:AppData/DiffStudio/library/pollution.ivp
//editor: Compute2:RichFunctionViewEditor
//input: double _t0 = 0 {units: min; caption: Initial; category: Time; min: 0; max: 0.9} [Initial time of simulation]
//input: double _t1 = 60 {units: min; caption: Final; category: Time; min: 1; max: 100; step: 1} [Final time of simulation]
//input: double _h = 0.1 {units: min; caption: Step; category: Time; min: 0.001; max: 0.1; step: 0.001} [Time step of simulation]
//input: double y1 = 0 {caption: NO2; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of NO2]
//input: double y2 = 0.2 {caption: NO; category: Initial concentrations; min: 0; max: 0.4} [Initial concentration of NO]
//input: double y3 = 0 {caption: O3P; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of O3P]
//input: double y4 = 0.04 {caption: O3; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of O3]
//input: double y5 = 0 {caption: HO2; category: Initial concentrations} [Initial concentration of HO2]
//input: double y6 = 0 {caption: OH; category: Initial concentrations} [Initial concentration of OH]
//input: double y7 = 0.1 {caption: HCHO; category: Initial concentrations} [Initial concentration of HCHO]
//input: double y8 = 0.3 {caption: CO; category: Initial concentrations} [Initial concentration of CO]
//input: double y9 = 0.01 {caption: ALD; category: Initial concentrations} [Initial concentration of ALD]
//input: double y10 = 0 {caption: MEO2; category: Initial concentrations} [Initial concentration of MEO2]
//input: double y11 = 0 {caption: C2O3; category: Initial concentrations} [Initial concentration of C2O3]
//input: double y12 = 0 {caption: CO2; category: Initial concentrations} [Initial concentration of CO2]
//input: double y13 = 0 {caption: PAN; category: Initial concentrations} [Initial concentration of PAN]
//input: double y14 = 0 {caption: CH3O; category: Initial concentrations} [Initial concentration of CH3O]
//input: double y15 = 0 {caption: HNO3; category: Initial concentrations} [Initial concentration of HNO3]
//input: double y16 = 0 {caption: O1D; category: Initial concentrations} [Initial concentration of O1D]
//input: double y17 = 0.007 {caption: SO2; category: Initial concentrations} [Initial concentration of SO2]
//input: double y18 = 0 {caption: SO4; category: Initial concentrations} [Initial concentration of SO4]
//input: double y19 = 0 {caption: NO3; category: Initial concentrations} [Initial concentration of NO3]
//input: double y20 = 0 {caption: N2O5; category: Initial concentrations} [Initial concentration of N2O5]
//input: double k1 = 0.35 {category: Reaction constants} [NO2 -> NO + O3P]
//input: double k2 = 26.6 {category: Reaction constants} [NO + O3 -> NO2]
//input: double k3 = 12300 {category: Reaction constants} [HO2 + NO -> NO2 + OH]
//input: double k4 = 0.00086 {category: Reaction constants} [HCHO -> 2 HO2 + CO]
//input: double k5 = 0.00082 {category: Reaction constants} [HCHO -> CO]
//input: double k6 = 15000 {category: Reaction constants} [HCHO + OH -> HO2 + CO]
//input: double k7 = 0.00013 {category: Reaction constants} [ALD -> MEO2 + HO2 + CO]
//input: double k8 = 24000 {category: Reaction constants} [ALD + OH -> C2O3]
//input: double k9 = 16500 {category: Reaction constants} [C2O3 + NO-> NO2 + MEO2 + CO2]
//input: double k10 = 9000 {category: Reaction constants} [C2O3 + NO2-> PAN]
//input: double k11 = 0.022 {category: Reaction constants} [PAN-> CH3O + NO2]
//input: double k12 = 12000 {category: Reaction constants} [MEO2 + NO-> CH3O + NO2]
//input: double k13 = 1.88 {category: Reaction constants} [CH3O-> HCHO + HO2]
//input: double k14 = 16300 {category: Reaction constants} [NO2 + OH -> HNO3]
//input: double k15 = 4800000 {category: Reaction constants} [O3P -> O3]
//input: double k16 = 0.00035 {category: Reaction constants} [O3 -> O1D]
//input: double k17 = 0.0175 {category: Reaction constants} [O3 -> O3P]
//input: double k18 = 100000000 {category: Reaction constants} [O1D -> 2 OH]
//input: double k19 = 444000000000 {category: Reaction constants} [O1D -> O3P]
//input: double k20 = 1240 {category: Reaction constants} [SO2 + OH -> SO4 + HO2]
//input: double k21 = 2.1 {category: Reaction constants} [NO3 -> NO]
//input: double k22 = 5.78 {category: Reaction constants} [NO3 -> NO2 + O3P]
//input: double k23 = 0.0474 {category: Reaction constants} [NO2 + O3 -> NO3]
//input: double k24 = 1780 {category: Reaction constants} [NO3 + NO2 -> N2O5]
//input: double k25 = 3.12 {category: Reaction constants} [N2O5 -> NO3 + NO2]
//output: dataframe df {caption: Pollution; viewer: Grid(block: 100) | Line chart(block: 100, multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | DiffStudio Facet(block: 100)}
export async function ivpModel_Pollution(_t0: number, _t1: number, _h: number, y1: number, y2: number, y3: number, y4: number, y5: number, y6: number, y7: number, y8: number, y9: number, y10: number, y11: number, y12: number, y13: number, y14: number, y15: number, y16: number, y17: number, y18: number, y19: number, y20: number, k1: number, k2: number, k3: number, k4: number, k5: number, k6: number, k7: number, k8: number, k9: number, k10: number, k11: number, k12: number, k13: number, k14: number, k15: number, k16: number, k17: number, k18: number, k19: number, k20: number, k21: number, k22: number, k23: number, k24: number, k25: number): Promise<DG.DataFrame> {
  return await runDiffStudioModel('System:AppData/DiffStudio/library/pollution.ivp', {_t0, _t1, _h, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25});
}

//name: Bioreactor
//description: Bioreactor simulation demo - multi-stage ODE with UF/DF mode switching. Case study: controlled Fab-arm exchange (cFAE) for bispecific antibody assembly
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/_bioreactor.png
//meta.help: models/bioreactor.md
//meta.role: model
//meta.diffStudioModel: System:AppData/DiffStudio/models/bioreactor.ivp
//editor: Compute2:RichFunctionViewEditor
//input: string mode {caption: Scenario; category: Process parameters; choices: OpenFile("System:AppData/DiffStudio/library/bioreactor-inputs.csv"); propagateChoice: all} [Load a scenario, then experiment — change MEA dose, temperature, switch time, or any other input]
//input: double _t0 = 0 {units: min; caption: Initial; category: Misc}                       [Simulation start time]
//input: double _t1 = 200 {units: min; caption: Reduction;   category: Duration; min: 200; max: 500}  [Duration of the cysteamine reduction stage (cFAE half-Ab exchange)]
//input: double _h = 0.5 {units: min; caption: Step;    category: Misc; min: 0.1; max: 2}     [ODE solver time step. Smaller = more accurate]
//input: double FFox = 0.2 {units: mmol/L; category: Initial values; min: 0.15; max: 0.25; step: 0.01}  [Parental F405L IgG1 (intact homodimer) — initial concentration]
//input: double KKox = 0.2 {units: mmol/L; category: Initial values; min: 0.15; max: 0.25; step: 0.01}  [Parental K409R IgG1 (intact homodimer) — initial concentration]
//input: double FFred = 0.1 {units: mmol/L; category: Initial values; min: 0.08; max: 0.12; step: 0.01}  [F405L homodimer with hinge disulfides reduced]
//input: double KKred = 0.1 {units: mmol/L; category: Initial values; min: 0.08; max: 0.12; step: 0.01}  [K409R homodimer with hinge disulfides reduced]
//input: double Ffree = 0 {units: mmol/L; category: Initial values}                                    [F405L half-antibody (exchange-competent intermediate)]
//input: double Kfree = 0 {units: mmol/L; category: Initial values}                                    [K409R half-antibody (exchange-competent intermediate)]
//input: double FKred = 0 {units: mmol/L; category: Initial values}                                    [Bispecific F/K heterodimer, hinge still reduced]
//input: double FKox = 0 {units: mmol/L; category: Initial values}                                    [Bispecific product with re-oxidized hinges (final cFAE product)]
//input: double MEAthiol = 15 {units: mmol/L; category: Initial values; min: 10;   max: 75}                [Cysteamine (2-MEA) dose. Typical cFAE: 25–75 mM]
//input: double DO2 = 0.12 {units: mmol/L; category: Initial values; min: 0.09; max: 2}                 [Dissolved O₂ concentration]
//input: double yO2P = 0.209 {units: atm;    category: Initial values}                                    [Headspace O₂ partial pressure. 0.209 ≈ ambient air]
//input: double CYST = 0 {units: mmol/L; category: Initial values}                                    [Cystamine — oxidized form of 2-MEA, builds up via air oxidation]
//input: double VL = 7.2 {units: L;      category: Initial values}                                    [Reactor liquid volume]
//input: double filtration = 300 {caption: Filtration; min: 100; max: 500; units: min; category: Duration} [Duration of UF/DF stage (concentration + diafiltration to remove 2-MEA)]
//input: double qin = 1 {units: L/min; caption: Gas;         category: Parameters;  min: 0.5; max: 2}              [Headspace gas sweep flow rate]
//input: double yO2in = 0.21 {              caption: O2 fraction; category: Parameters;  min: 0.1; max: 0.9}            [O₂ mole fraction in inlet gas. 0.21 = air, 1.0 = pure O₂]
//input: double T = 300 {units: K;     caption: Temperature; category: Parameters;  min: 250; max: 350}            [Reactor temperature (K). 300 K ≈ 27 °C]
//input: double P = 1 {units: atm;   caption: Pressure;    category: Parameters;  min: 1;   max: 2}              [Total headspace pressure]
//input: double switchTime = 135 {units: min;   caption: Switch;   category: Duration;        min: 70;  max: 180; step: 10}  [Transition from UF concentration (VL falls) to diafiltration (VL constant)]
//output: dataframe df {caption: Bioreactor; viewer: Grid(block: 100) | Line chart(block: 100, multiAxis: "true", segmentColumnName: "_Stage", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | DiffStudio Facet(block: 100, segmentColumnName: "_Stage")}
export async function ivpModel_Bioreactor(mode: string, _t0: number, _t1: number, _h: number, FFox: number, KKox: number, FFred: number, KKred: number, Ffree: number, Kfree: number, FKred: number, FKox: number, MEAthiol: number, DO2: number, yO2P: number, CYST: number, VL: number, filtration: number, qin: number, yO2in: number, T: number, P: number, switchTime: number): Promise<DG.DataFrame> {
  return await runDiffStudioModel('System:AppData/DiffStudio/models/bioreactor.ivp', {_t0, _t1, _h, FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, FKox, MEAthiol, DO2, yO2P, CYST, VL, filtration, qin, yO2in, T, P, switchTime});
}

//name: PK-PD
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/pkpd.png
//meta.help: models/pk-pd.md
//meta.role: model
//meta.diffStudioModel: System:AppData/DiffStudio/models/pk-pd.ivp
//editor: Compute2:RichFunctionViewEditor
//input: int _count = 10 {caption: Count; category: Dosing; min: 1; max: 20} [Number of doses delivered in the loop]
//input: double _t0 = 0 {units: h; caption: Begin;    category: Misc;   min: 0; max: 1}     [Start of the dosing interval]
//input: double _t1 = 12 {units: h; caption: Interval; category: Dosing; min: 5; max: 15}    [Length of one dosing interval (hours between doses)]
//input: double _h = 0.2 {units: h; caption: Step;    category: Misc;   min: 0.1; max: 2}   [ODE solver time step]
//input: double depot = 0 {caption: Depo; category: Misc}                     [Amount in the absorption site at t=0]
//input: double centr = 0 {caption: Central;     category: Misc}              [Amount in the central compartment at t=0]
//input: double peri = 0 {caption: Peripheral;  category: Misc}              [Amount in the peripheral compartment at t=0]
//input: double eff = 0.2 {caption: Init effect; category: Misc}              [Initial effect / biomarker level]
//input: double C2 = 0 {caption: Central;     category: Initial concentrations}  [Central concentration at t=0]
//input: double C3 = 0 {caption: Peripheral;  category: Initial concentrations}  [Peripheral concentration at t=0]
//input: double dose = 10000 {caption: Dose; category: Dosing; min: 1e3; max: 2e4; step: 1e3}             [Amount given at each dose]
//input: double KA = 0.3 {caption: Rate constant; category: PK parameters; min: 0.1; max: 1}          [Absorption rate constant (1/h)]
//input: double CL = 2 {caption: Clearance;     category: PK parameters; min: 1;   max: 5}          [Total clearance from central compartment]
//input: double V2 = 4 {caption: Central volume; category: PK parameters; min: 1;   max: 10}        [Volume of the central compartment]
//input: double Q = 1 {caption: Inter rate;    category: PK parameters; min: 0.1; max: 1}          [Inter-compartmental clearance between central and peripheral]
//input: double V3 = 30 {caption: Peri volume;   category: PK parameters; min: 20;  max: 40}         [Volume of the peripheral compartment]
//input: double EC50 = 8 {caption: Effect;        category: PD parameters; min: 1;   max: 10}         [Drug concentration giving half-maximal effect]
//input: double Rate = 0.2 {category: PD parameters; min: 0.1; max: 0.5}                                [Turnover rate of the effect variable]
//output: dataframe df {caption: PK-PD; viewer: Grid(block: 100) | Line chart(block: 100, multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | DiffStudio Facet(block: 100)}
export async function ivpModel_PK_PD(_count: number, _t0: number, _t1: number, _h: number, depot: number, centr: number, peri: number, eff: number, C2: number, C3: number, dose: number, KA: number, CL: number, V2: number, Q: number, V3: number, EC50: number, Rate: number): Promise<DG.DataFrame> {
  return await runDiffStudioModel('System:AppData/DiffStudio/models/pk-pd.ivp', {_count, _t0, _t1, _h, depot, centr, peri, eff, C2, C3, dose, KA, CL, V2, Q, V3, EC50, Rate});
}
