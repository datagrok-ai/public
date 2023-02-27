/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

import {callWasm} from '../wasm/callWasm';
import {getCppInput, getResult} from '../wasm/callWasmUtils';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initFAEexplicit();
}

//name: solveFAEexplicit
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 1000.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _CystamineInitial = 0.0 {units: ; caption: Cystamine; category: initial values}
//input: double _FFoxInitial = 0.268 {units: ; caption: FFox; category: initial values}
//input: double _FFredInitial = 0.0 {units: ; caption: FFred; category: initial values}
//input: double _KKredInitial = 0.0 {units: ; caption: KKred; category: initial values}
//input: double _KfreeInitial = 0.0 {units: ; caption: Kfree; category: initial values}
//input: double _FKredInitial = 0.0 {units: ; caption: FKred; category: initial values}
//input: double _KKoxInitial = 0.268 {units: ; caption: KKox; category: initial values}
//input: double _FKoxInitial = 0.0 {units: ; caption: FKox; category: initial values}
//input: double _MEAthiolInitial = 34.0 {units: ; caption: MEAthiol; category: initial values}
//input: double _CO2Initial = 0.22 {units: ; caption: CO2; category: initial values}
//input: double _yO2PInitial = 0.209 {units: ; caption: yO2P; category: initial values}
//input: double _VLInitial = 6.6 {units: ; caption: VL; category: initial values}
//input: double _FfreeInitial = 0.0 {units: ; caption: Ffree; category: initial values}
//input: double _qinVal = 1.0 {units: Liters gas/minute; caption: qin; category: parameters}
//input: double _percentO2saturationVal = 100.0 {units: ; caption: percentO2saturation; category: parameters}
//input: double _yO2inVal = 0.209 {units: ; caption: yO2in; category: parameters}
//input: double _pKa2MEAVal = 8.19 {units: ; caption: pKa2MEA; category: parameters}
//input: double _HVal = 1.072069378 {units: ; caption: H; category: parameters}
//input: double _TVal = 300.0 {units: degK; caption: T; category: parameters}
//input: double _RVal = 0.082 {units: Liter Atm / mole degK; caption: R; category: parameters}
//input: double _PVal = 1.0 {units: atma; caption: P; category: parameters}
//output: dataframe dfSolution {caption: FAE(explicit); viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function solveFAEexplicit(initial, final, step,
  _CystamineInitial, _FFoxInitial, _FFredInitial, _KKredInitial, _KfreeInitial, _FKredInitial, _KKoxInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _VLInitial, _FfreeInitial, 
  _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal)
{
  let _tCount = Math.trunc((final - initial) / step) + 1;
  let _varsCount = 14;

  return callWasm(FAEexplicit, 'solveFAEexplicit',
    [ initial, final, step,
     _CystamineInitial, _FFoxInitial, _FFredInitial, _KKredInitial, _KfreeInitial, _FKredInitial, _KKoxInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _VLInitial, _FfreeInitial,
     _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal,
     _tCount, _varsCount ] );
}

//name: solveFAEexplicitWebWorker
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 1000.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _CystamineInitial = 0.0 {units: ; caption: Cystamine; category: initial values}
//input: double _FFoxInitial = 0.268 {units: ; caption: FFox; category: initial values}
//input: double _FFredInitial = 0.0 {units: ; caption: FFred; category: initial values}
//input: double _KKredInitial = 0.0 {units: ; caption: KKred; category: initial values}
//input: double _KfreeInitial = 0.0 {units: ; caption: Kfree; category: initial values}
//input: double _FKredInitial = 0.0 {units: ; caption: FKred; category: initial values}
//input: double _KKoxInitial = 0.268 {units: ; caption: KKox; category: initial values}
//input: double _FKoxInitial = 0.0 {units: ; caption: FKox; category: initial values}
//input: double _MEAthiolInitial = 34.0 {units: ; caption: MEAthiol; category: initial values}
//input: double _CO2Initial = 0.22 {units: ; caption: CO2; category: initial values}
//input: double _yO2PInitial = 0.209 {units: ; caption: yO2P; category: initial values}
//input: double _VLInitial = 6.6 {units: ; caption: VL; category: initial values}
//input: double _FfreeInitial = 0.0 {units: ; caption: Ffree; category: initial values}
//input: double _qinVal = 1.0 {units: Liters gas/minute; caption: qin; category: parameters}
//input: double _percentO2saturationVal = 100.0 {units: ; caption: percentO2saturation; category: parameters}
//input: double _yO2inVal = 0.209 {units: ; caption: yO2in; category: parameters}
//input: double _pKa2MEAVal = 8.19 {units: ; caption: pKa2MEA; category: parameters}
//input: double _HVal = 1.072069378 {units: ; caption: H; category: parameters}
//input: double _TVal = 300.0 {units: degK; caption: T; category: parameters}
//input: double _RVal = 0.082 {units: Liter Atm / mole degK; caption: R; category: parameters}
//input: double _PVal = 1.0 {units: atma; caption: P; category: parameters}
export async function solveFAEexplicitWebWorker(initial, final, step,
  _CystamineInitial, _FFoxInitial, _FFredInitial, _KKredInitial, _KfreeInitial, _FKredInitial, _KKoxInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _VLInitial, _FfreeInitial, 
  _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal)
{
  let _tCount = Math.trunc((final - initial) / step) + 1;
  let _varsCount = 14;

  // web worker declaration
  var worker = new Worker(new URL('../wasm/workerFAE.js', import.meta.url));              

  // post function specification & inputs
  worker.postMessage(getCppInput(FAEexplicit['solveFAEexplicit'].arguments, 
    [ initial, final, step,
      _CystamineInitial, _FFoxInitial, _FFredInitial, _KKredInitial, _KfreeInitial, _FKredInitial, _KKoxInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _VLInitial, _FfreeInitial,
      _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal,
      _tCount, _varsCount ]));       
  
  worker.onmessage = function(e) {       
    // get results of wasm-computations 
    let output = getResult(FAEexplicit['solveFAEexplicit'], e.data);
  
    output.name = 'FAE(explicit, webworker)';
    let view = grok.shell.addTableView(output);
    view.lineChart({ markerType: 'dot', sharex: 'true', multiAxis: 'true'});
  }    
}