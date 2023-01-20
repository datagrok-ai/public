//name: Default simple script JS
//description: Shows that default view works (or not)
//language: javascript
//tags: bioreactor
//input: double S1 = 1 {caption: S1; units: ft² * BTU / hr °F; category: Category1}
//input: double S2 = 100 {caption: Initial scalar; units: L/min; category: Category1}
//input: double S3 = -2 {caption: Dummy scalar.; units: °C; category: Category1}
//input: double S4 = 5 {caption: Final S4; units: cells/mL; category: Category2}
//input: double S5 = 210.5 {caption: Final Volume; units: L; category: Category2}
//input: double S = 31.1 {caption: Temp.; units: °C; category: Category2}
//output: dataframe tempOnTime { viewer: Line chart(block: 50, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Statistics(block: 50); category: CoolingRate }
//output: dataframe tempOnTime2 { viewer: Line chart(block: 25, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Grid(block: 75); category: CoolingRate2 }
//output: double O1 {caption: Temp 1.; units: °C; category: CoolingRate2 }
//output: dataframe tempOnTime3 { caption: "My lovely caption"; viewer: Line chart(block: 75, x: "Time (hours)", y: "Temperature (°C)", showSplitSelector: false) | Grid(block: 25); category: CoolingRate2 }
//output: double O2 {caption: Temp 2.; units: °C; category: CoolingRate3 }
//output: double O3 {caption: Temp 3.; units: °C; category: CoolingRate3 }
//output: double O4 {caption: Temp 4.; units: °C; category: CoolingRate3 }
//output: double O5 {caption: Temp 5.; units: °C; category: CoolingRate3 }
//editor: Compute:RichFunctionViewEditor

const tt = [...Array(S2).keys()]
const Sol = [...Array(S2).keys()].map((x) => x*2)

tempOnTime = DG.DataFrame.fromColumns([
    DG.Column.float("Time (hours)", S2).init((i) => tt[i]),
    DG.Column.float("Temperature (°C)", S2).init((i) => Sol[i]),
])
tempOnTime2 = tempOnTime
tempOnTime3 = tempOnTime

O1 = S2
O2 = S2
O3 = S2
O4 = S2
O5 = S2