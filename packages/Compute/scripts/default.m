#name: Default simple script
#description: Shows that default view works (or not)
#language: octave
#tags: bioreactor
#input: double S1 = 424.2832 {caption: S1; units: ft² * BTU / hr °F}
#input: double S2 = 12 {caption: Initial scalar; units: L/min}
#input: double S3 = -2 {caption: Dummy scalar.; units: °C}
#input: double S4 = 5 {caption: Final S4; units: cells/mL}
#input: double S5 = 210.5 {caption: Final Volume; units: L}
#input: double S = 31.1 {caption: Temp.; units: °C}
#output: dataframe tempOnTime { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate }
#editor: Compute:FunctionViewEditor

tt = zeros(10, 1)
Sol = ones(10, 1)

labels = {"Time (hours)", "Temperature (°C)"};
tempOnTime = [labels; num2cell([tt, Sol])];