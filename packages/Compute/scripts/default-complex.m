#name: Default complex script
#description: Shows that default view works (or not)
#language: octave
#tags: bioreactor
#input: dataframe test { viewer: Line chart(x: "Axis 1", y: "Axis2", showSplitSelector: false) }
#input: double S1 = 424.2832 {caption: S1; units: ft² * BTU / hr °F}
#input: double S2 = 12 {caption: Initial scalar; units: L/min}
#input: double S3 = -2 {caption: Dummy scalar.; units: °C}
#input: double S4 = 5 {caption: Final S4; units: cells/mL}
#input: double S5 = 210.5 {caption: Final Volume; units: L}
#input: double S = 31.1 {caption: Temp.; units: °C}
#output: dataframe tempOnTime { viewer: Line chart(x: "Axis 1", y: "Axis2", showSplitSelector: false); category: Cooling rate }
#output: double IC2 {caption: Culture temp.; units: °C; category: Cooling rate}
#output: dataframe tempOnTime2 { viewer: Line chart(x: "Axis 1", y: "Axis2", showSplitSelector: false) | Grid(block: 50); category: Cooling rate 2 }
#output: dataframe tempOnTime3 { viewer: Line chart(block: 50, x: "Axis 1", y: "Axis2", showSplitSelector: false) | Line chart(block: 50, x: "Axis 1", y: "Axis2", showSplitSelector: false); category: Cooling rate 3 }
#output: dataframe tempOnTime4 { viewer: Line chart(block: 50, x: "Axis 1", y: "Axis2", showSplitSelector: false); category: Cooling rate 3 }
#output: dataframe tempOnTime5 { viewer: Line chart(block: 50, x: "Axis 1", y: "Axis2", showSplitSelector: false) | Line chart(block: 50, x: "Axis 1", y: "Axis2", showSplitSelector: false); category: Cooling rate 3 }
#output: double IC3 {caption: Culture temp.; units: °C; category: Cooling rate 3}
#editor: Compute:FunctionViewEditor

tempOnTime = test;
tempOnTime2 = tempOnTime;
tempOnTime3 = tempOnTime;
tempOnTime4 = tempOnTime;
tempOnTime5 = tempOnTime;

IC2 = S1
IC3 = S1