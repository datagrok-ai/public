 // THIS FILE IS GENERATED AUTOMATICALLY. DO NOT CHANGE ANYTHING!

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

import { callWasm } from '../wasm/callWasm';

//tags: init
export async function init() {
  await initFAEextended();
}

//name: solveFAEextended
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 1000.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _FFoxInitial = 0.268 {units: ; caption: FFox; category: initial values}
//input: double _KKoxInitial = 0.268 {units: ; caption: KKox; category: initial values}
//input: double _FFredInitial = 0.0 {units: ; caption: FFred; category: initial values}
//input: double _KKredInitial = 0.0 {units: ; caption: KKred; category: initial values}
//input: double _FfreeInitial = 0.0 {units: ; caption: Ffree; category: initial values}
//input: double _KfreeInitial = 0.0 {units: ; caption: Kfree; category: initial values}
//input: double _FKredInitial = 0.0 {units: ; caption: FKred; category: initial values}
//input: double _FKoxInitial = 0.0 {units: ; caption: FKox; category: initial values}
//input: double _MEAthiolInitial = 34.0 {units: ; caption: MEAthiol; category: initial values}
//input: double _CO2Initial = 0.22 {units: ; caption: CO2; category: initial values}
//input: double _yO2PInitial = 0.209 {units: ; caption: yO2P; category: initial values}
//input: double _CystamineInitial = 0.0 {units: ; caption: Cystamine; category: initial values}
//input: double _VLInitial = 6.6 {units: ; caption: VL; category: initial values}
//input: double _qinVal = 1.0 {units: Liters gas/minute; caption: qin; category: parameters}
//input: double _percentO2saturationVal = 100.0 {units: ; caption: percentO2saturation; category: parameters}
//input: double _yO2inVal = 0.209 {units: ; caption: yO2in; category: parameters}
//input: double _pKa2MEAVal = 8.19 {units: ; caption: pKa2MEA; category: parameters}
//input: double _HVal = 1.072069378 {units: ; caption: H; category: parameters}
//input: double _TVal = 300.0 {units: degK; caption: T; category: parameters}
//input: double _RVal = 0.082 {units: Liter Atm / mole degK; caption: R; category: parameters}
//input: double _PVal = 1.0 {units: atma; caption: P; category: parameters}
//input: double _TimeToSwitchVal = 180.0 {units: ; caption: TimeToSwitch; category: parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function solveFAEextended(initial, final, step,
  _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _CystamineInitial, _VLInitial, 
  _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal)
{
  let _tCount = Math.trunc((final - initial) / step) + 1;
  let _varsCount = 14;

  return callWasm(FAEextended, 'solveFAEextended',
    [ initial, final, step,
     _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _CystamineInitial, _VLInitial,
     _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal,
     _tCount, _varsCount ] );
}
