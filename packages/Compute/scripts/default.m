#name: Default simple script
#description: Shows that default view works (or not)
#language: octave
#tags: bioreactor
#input: double UA = 424.2832 {caption: UA; units: ft² * BTU / hr °F}
#input: double CoolingWaterFlow = 12 {caption: Coolant flowrate; units: L/min}
#input: double Tcoolingin = -2 {caption: Coolant temp.; units: °C}
#input: double VCD = 5 {caption: Target VCD; units: cells/mL}
#input: double VL = 210.5 {caption: Final Volume; units: L}
#input: double IC = 31.1 {caption: Culture temp.; units: °C}
#output: dataframe tempOnTime { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate }
#editor: Compute:FunctionViewEditor

pkg load symbolic

## Conversion
CalperBTU = 252;
CalperhourperWatt = 3600/4.183076;
Cp = 1;
CoolingFactor = 56.5;

## Explicit Equation
##specificenthalpyreaction = 3.7136 * exp(0.0652*Tsub);
##Qrxn = VCD * 10^6 * VL * specificenthalpyofreaction/10^9;
##LMTD = (Tcoolingin-(Tsub + (Tcoolingin-Tsub)/(exp(UA*(CoolingFactor/100)*1.8*CalperBTU/(CoolingWaterFlow*Cp*60*1000)))))/log((Tsub-(Tsub + (Tcoolingin-Tsub)/(exp(UA*(CoolingFactor/100)*1.8*CalperBTU/(CoolingWaterFlow*Cp*60*1000)))))/(Tsub-Tcoolingin));
##Tcoolingout =  Tsub + (Tcoolingin-Tsub)/(exp(UA*(CoolingFactor/100)*1.8*CalperBTU/(CoolingWaterFlow*Cp*60*1000)));
## Differential equation
syms Tsub;

F1 = @(t, Tsub) -UA * (CoolingFactor / 100) * ((Tcoolingin - (Tsub + (Tcoolingin - Tsub) / (exp(UA * (CoolingFactor / 100) * 1.8 * CalperBTU / (CoolingWaterFlow * Cp * 60 * 1000))))) / log((Tsub - (Tsub + (Tcoolingin - Tsub) / (exp(UA * (CoolingFactor / 100) * 1.8 * CalperBTU / (CoolingWaterFlow * Cp * 60 * 1000))))) / (Tsub - Tcoolingin))) * 1.8 * CalperBTU / (VL * Cp * 1000) + (VCD * 10^6 * VL * 3.7136 * exp(0.0652 * Tsub) / 10^9) * CalperhourperWatt / (VL * Cp * 1000);

## Change the time here
t = 3;
dt = 3/20;

## Solving the differential equations
[tt, Sol] = ode45(F1, 0:dt:t, IC);
labels = {"Time (hours)", "Temperature (°C)"};
tempOnTime = [labels; num2cell([tt, Sol])];