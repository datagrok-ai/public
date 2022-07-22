#name: Default simple script
#description: Shows that default view works (or not)
#language: octave
#tags: bioreactor
#input: double AA = 424.2832 {caption: AA; units: ft² * BTU / hr °F}
#input: double FWC = 12 {caption:FWC; units: L/min}
#input: double Tcoolingin = -2 {caption: Coolant temp.; units: °C}
#input: double VCD = 5 {caption: Target VCD; units: cells/mL}
#input: double VL = 210.5 {caption: Final Volume; units: L}
#input: double IC = 31.1 {caption: Temp.; units: °C}
#output: dataframe tempOnTime { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate }
#editor: Compute:FunctionViewEditor

tt = zeros(10, 1)
Sol = ones(10, 1)

labels = {"Time (hours)", "Temperature (°C)"};
tempOnTime = [labels; num2cell([tt, Sol])];