#name: Default complex script
#description: Shows that default view works (or not)
#language: octave
#tags: bioreactor
#input: dataframe test { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) }
#input: double AA = 424.2832 {caption: AA; units: ft² * BTU / hr °F}
#input: double FWC = 12 {caption: Bblabla flowrate; units: L/min}
#input: double Dummy scalar = -2 {caption: Dummy scalar.; units: °C}
#input: double VCD = 5 {caption: VCD; units: cells/mL}
#input: double VL = 210.5 {caption: Final Volume; units: L}
#input: double IC = 31.1 {caption: Temp.; units: °C}
#output: dataframe tempOnTime { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate }
#output: double IC2 {caption: Culture temp.; units: °C; category: Cooling rate}
#output: dataframe tempOnTime2 { viewer: Line chart(x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Grid(block: 50); category: Cooling rate 2 }
#output: dataframe tempOnTime3 { viewer: Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate 3 }
#output: dataframe tempOnTime4 { viewer: Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate 3 }
#output: dataframe tempOnTime5 { viewer: Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false); category: Cooling rate 3 }
#output: double IC3 {caption: Culture temp.; units: °C; category: Cooling rate 3}
#editor: Compute:FunctionViewEditor

tempOnTime = test;
tempOnTime2 = tempOnTime;
tempOnTime3 = tempOnTime;
tempOnTime4 = tempOnTime;
tempOnTime5 = tempOnTime;

IC2 = IC
IC3 = IC