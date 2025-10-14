//name: Object cooling
//description: Uses Newton's law of cooling to simulate object cooling process. Default values are for cube of boiling water in air.
//language: javascript
//tags: simulation, demo
//input: double ambTemp = 22 {caption: Ambient temperature; units: C; category: Environment; block: 50; validatorFunc: Compute:AmbTempValidator; }
//input: double initTemp = 100 {caption: Initial temperature; units: C; category: Environment; block: 50; validatorFunc: Compute:InitialTempValidator; }
//input: double desiredTemp = 30 {caption: Desired temperature; units: C; category: Environment; validatorFunc: Compute:DesiredTempValidator; }
//input: double area = 0.06 {caption: Surface area; units: m²; category: Object properties}
//input: double heatCap = 4200 {caption: Heat capacity; units: J/C; category: Object properties; validatorFunc: Compute:HeatCapValidator; }
//input: double heatTransferCoeff = 8.3 {caption: Heat transfer coefficient; units: W/(m² * C); category: Object properties}
//input: string previousRun {caption: Previous run; category: Simulation; input: Compute:ObjectCoolingSelector; funccallMapping: {"simTime": "timeToCool"}; nullable: true}
//input: int simTime = 21600 {caption: Simulation time; units: sec; category: Simulation; validatorFunc: Compute:SimTimeValidator; validatorFuncOptions: { "reasonableMin": 10800, "reasonableMax": 100000} }
//output: dataframe simulation {caption: Temp. vs time; category: Output; viewer: Line chart(block: 75) | Grid(block: 25)}
//output: double timeToCool {caption: Time to cool; units: sec.; category: Output}
//output: double coolingFactor {caption: Cooling factor; units: 1 / sec.; category: Calculations; precision: 5}
//output: double tempDiff {caption: Temperature difference; units: C; category: Calculations}
//editor: Compute:RichFunctionViewEditor
//meta.features: {"sens-analysis": true, "upload": true}
//meta.foldedCategories: ["Object properties"]
//meta.compareCustomizer: Compute:CustomCustomizer
//meta.uploadFunc: Compute:CustomUploader
//meta.help: Compute/readme.md
//meta.mainParams: ["ambTemp", "initTemp", "desiredTemp"]

timeToCool = undefined;

const tempDiff = initTemp - ambTemp;
const coolingFactor = heatTransferCoeff * area / heatCap;

const timeStamps = new Float32Array(simTime).map((_, idx) => idx);
const simulatedTemp = timeStamps.map((timeStamp) => {
  const currentTemp = ambTemp + (tempDiff * (Math.E ** -(coolingFactor * timeStamp)));

  if (!timeToCool && currentTemp < desiredTemp) {
    timeToCool = timeStamp;
  }

  return currentTemp;
});

simulation = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('Time', timeStamps),
  DG.Column.fromFloat32Array('Temperature', simulatedTemp),
]);
