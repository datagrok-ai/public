/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {_initReactor, _simulateBioreactor} from '../wasm/reactorAPI';
import {customRun, showCustomRunResults} from '../wasm/demoTools';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() { 
  await _initReactor(); 
}

//name: Bioreactor
//description: Bioreactor simulation.
//tags: model
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 1000.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _FFoxInitial = 0.268 {units: mMole/Liter; caption: FFox; category: initial values}
//input: double _KKoxInitial = 0.268 {units: mMole/Liter; caption: KKox; category: initial values}
//input: double _FFredInitial = 0.0 {units: mMole/Liter; caption: FFred; category: initial values}
//input: double _KKredInitial = 0.0 {units: mMole/Liter; caption: KKred; category: initial values}
//input: double _FfreeInitial = 0.0 {units: mMole/Liter; caption: Ffree; category: initial values}
//input: double _KfreeInitial = 0.0 {units: mMole/Liter; caption: Kfree; category: initial values}
//input: double _FKredInitial = 0.0 {units: mMole/Liter; caption: FKred; category: initial values}
//input: double _FKoxInitial = 0.0 {units: mMole/Liter; caption: FKox; category: initial values}
//input: double _MEAthiolInitial = 34.0 {units: mMole/Liter; caption: MEAthiol; category: initial values}
//input: double _CO2Initial = 0.22 {units: ; caption: CO2; category: initial values}
//input: double _yO2PInitial = 0.209 {units: ATMa O2; caption: yO2P; category: initial values}
//input: double _CystamineInitial = 0.0 {units: mMole/Liter; caption: Cystamine; category: initial values}
//input: double _VLInitial = 6.6 {units: Liters; caption: VL; category: initial values}
//input: double _qinVal = 1.0 {units: Liters gas/minute; caption: qin; category: parameters}
//input: double _percentO2saturationVal = 100.0 {units: Dissolved Oxygen % saturation; caption: percentO2saturation; category: parameters}
//input: double _yO2inVal = 0.209 {units: mole fraction O2; caption: yO2in; category: parameters}
//input: double _pKa2MEAVal = 8.19 {units: pKa of MAB thiol to thiolate; caption: pKa2MEA; category: parameters}
//input: double _HVal = 1.072069378 {units: mmole O2/ L liquid /ATM O2; caption: H; category: parameters}
//input: double _TVal = 300.0 {units: degK; caption: T; category: parameters}
//input: double _RVal = 0.082 {units: Liter Atm / mole degK; caption: R; category: parameters}
//input: double _PVal = 1.0 {units: atma; caption: P; category: parameters}
//input: double _TimeToSwitchVal = 180.0 {units: minute; caption: TimeToSwitch; category: parameters}
//output: dataframe solution {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function Bioreactor(initial: number, final: number, step: number,
  _FFoxInitial: number, _KKoxInitial: number, _FFredInitial: number, _KKredInitial: number, 
  _FfreeInitial: number, _KfreeInitial: number, _FKredInitial: number, _FKoxInitial: number,
  _MEAthiolInitial: number, _CO2Initial: number, _yO2PInitial: number, _CystamineInitial: number, 
  _VLInitial: number, _qinVal: number, _percentO2saturationVal: number, _yO2inVal: number, 
  _pKa2MEAVal: number, _HVal: number, _TVal: number, _RVal: number, _PVal: number, 
  _TimeToSwitchVal: number): Promise<DG.DataFrame>
{
  return await _simulateBioreactor(initial, final, step,
    _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, 
    _FfreeInitial, _KfreeInitial, _FKredInitial, _FKoxInitial,
    _MEAthiolInitial, _CO2Initial, _yO2PInitial, _CystamineInitial, 
    _VLInitial, _qinVal, _percentO2saturationVal, _yO2inVal, 
    _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal, 
    _TimeToSwitchVal);
}

//name: Bioreactor Demo
//description: Bioreactor simulation demo.
//meta.demoPath: Bioreactors | Simulate bioreactor
export async function demoScript(): Promise<any>  {
  const demoScript = new DemoScript('Bioreactor Demo', 'Bioreactor simulation demo.'); 

  let tables: any;
    
  await demoScript    
    .step('Run', async () => { tables = await customRun() }, {description: 'Launch bioreactor simulation with different parameters.', delay: 5000})  
    .step('Compare', async () => { showCustomRunResults(tables) }, {description: 'Compare results of different runs.'})  
    .start();
}
