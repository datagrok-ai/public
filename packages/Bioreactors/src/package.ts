/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {_initinitBioreactor, _simulateBioreactor} from '../wasm/BioreactorAPI';
import {customRun, showCustomRunResults} from '../wasm/demoTools';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() { 
  await _initinitBioreactor(); 
}

//name: Bioreactor
//description: Controlled fab-arm exchange mechanism simulation.
//tags: model
//input: double initial = 0.0 {caption: Initial; category: Time, min}
//input: double final = 1000.0 {caption: Final; category: Time, min}
//input: double step = 0.1 {caption: Step; category: Time, min}
//input: double _FFoxInitial = 0.2 {units: mmol/L; caption: FF oxidized (FFox); category: Initial values} 
//input: double _KKoxInitial = 0.2 {units: mmol/L; caption: KK oxidized (KKox); category: Initial values}
//input: double _FFredInitial = 0.1 {units: mmol/L; caption: FF reduced (FFred); category: Initial values}
//input: double _KKredInitial = 0.1 {units: mmol/L; caption: KK reduced (KKred); category: Initial values}
//input: double _FfreeInitial = 0.0 {units: mmol/L; caption: F free (Ffree); category: Initial values}
//input: double _KfreeInitial = 0.0 {units: mmol/L; caption: K free (Kfree); category: Initial values}
//input: double _FKredInitial = 0.0 {units: mmol/L; caption: FK reduced (FKred); category: Initial values}
//input: double _FKoxInitial = 0.0 {units: mmol/L; caption: FK oxidized (FKox); category: Initial values}
//input: double _MEAthiolInitial = 15.0 {units: mmol/L; caption: MEAthiol (MEA); category: Initial values}
//input: double _CO2Initial = 0.12 {units: mmol/L; caption: Dissolved oxygen (CO2); category: Initial values}
//input: double _yO2PInitial = 0.209 {units: atm; caption: Atm headspace (yO2P); category: Initial values}
//input: double _CYSTInitial = 0.0 {units: mmol/L; caption: Cystamine (CYST); category: Initial values}
//input: double _VLInitial = 7.2 {units: L; caption: Liquid volume (VL); category: Initial values}
//input: double _qinVal = 1.0 {units: L/min; caption: Gas to headspace; category: Parameters}
//input: double _yO2inVal = 0.21 {units: ; caption: Oxygen mole fraction; category: Parameters}
//input: double _HVal = 1.3 {units: mmol/(L atm); caption: Henry's law constant; category: Parameters}
//input: double _TVal = 300.0 {units: K; caption: System temperature; category: Parameters}
//input: double _RVal = 0.082 {units: L atm/(mol K); caption: Gas constant; category: Parameters}
//input: double _PVal = 1.0 {units: atm; caption: Headspace pressure; category: Parameters}
//input: double _TimeToSwitchVal = 135.0 {units: min; caption: Switch mode time; category: Parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(block: 100, x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function Bioreactor(initial: number, final: number, step: number,
  _FFoxInitial: number, _KKoxInitial: number, _FFredInitial: number, _KKredInitial: number, 
  _FfreeInitial: number, _KfreeInitial: number, _FKredInitial: number, _FKoxInitial: number,
  _MEAthiolInitial: number, _CO2Initial: number, _yO2PInitial: number, _CYSTInitial: number, 
  _VLInitial: number, _qinVal: number, _yO2inVal: number, _HVal: number, _TVal: number, 
  _RVal: number, _PVal: number, _TimeToSwitchVal: number): Promise<DG.DataFrame>
{
  return await _simulateBioreactor(initial, final, step,
    _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, 
    _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, 
    _yO2PInitial, _CYSTInitial, _VLInitial, _qinVal, _yO2inVal, 
    _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal);
}

//name: Bioreactor Demo
//description: Controlled fab-arm exchange mechanism simulation.
//input: double initial = 0.0 {caption: Initial; category: Time, min}
//input: double final = 1000.0 {caption: Final; category: Time, min}
//input: double step = 0.1 {caption: Step; category: Time, min}
//input: double _FFoxInitial = 0.2 {units: mmol/L; caption: FF oxidized (FFox); category: Initial values} 
//input: double _KKoxInitial = 0.2 {units: mmol/L; caption: KK oxidized (KKox); category: Initial values}
//input: double _FFredInitial = 0.1 {units: mmol/L; caption: FF reduced (FFred); category: Initial values}
//input: double _KKredInitial = 0.1 {units: mmol/L; caption: KK reduced (KKred); category: Initial values}
//input: double _FfreeInitial = 0.0 {units: mmol/L; caption: F free (Ffree); category: Initial values}
//input: double _KfreeInitial = 0.0 {units: mmol/L; caption: K free (Kfree); category: Initial values}
//input: double _FKredInitial = 0.0 {units: mmol/L; caption: FK reduced (FKred); category: Initial values}
//input: double _FKoxInitial = 0.0 {units: mmol/L; caption: FK oxidized (FKox); category: Initial values}
//input: double _MEAthiolInitial = 15.0 {units: mmol/L; caption: MEAthiol (MEA); category: Initial values}
//input: double _CO2Initial = 0.12 {units: mmol/L; caption: Dissolved oxygen (CO2); category: Initial values}
//input: double _yO2PInitial = 0.209 {units: atm; caption: Atm headspace (yO2P); category: Initial values}
//input: double _CYSTInitial = 0.0 {units: mmol/L; caption: Cystamine (CYST); category: Initial values}
//input: double _VLInitial = 7.2 {units: L; caption: Liquid volume (VL); category: Initial values}
//input: double _qinVal = 1.0 {units: L/min; caption: Gas to headspace; category: Parameters}
//input: double _yO2inVal = 0.21 {units: ; caption: Oxygen mole fraction; category: Parameters}
//input: double _HVal = 1.3 {units: mmol/(L atm); caption: Henry's law constant; category: Parameters}
//input: double _TVal = 300.0 {units: K; caption: System temperature; category: Parameters}
//input: double _RVal = 0.082 {units: L atm/(mol K); caption: Gas constant; category: Parameters}
//input: double _PVal = 1.0 {units: atm; caption: Headspace pressure; category: Parameters}
//input: double _TimeToSwitchVal = 135.0 {units: min; caption: Switch mode time; category: Parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(block: 100, x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
//meta.runOnStart: true
export async function BioreactorDemo(initial: number, final: number, step: number,
  _FFoxInitial: number, _KKoxInitial: number, _FFredInitial: number, _KKredInitial: number, 
  _FfreeInitial: number, _KfreeInitial: number, _FKredInitial: number, _FKoxInitial: number,
  _MEAthiolInitial: number, _CO2Initial: number, _yO2PInitial: number, _CYSTInitial: number, 
  _VLInitial: number, _qinVal: number, _yO2inVal: number, _HVal: number, _TVal: number, 
  _RVal: number, _PVal: number, _TimeToSwitchVal: number): Promise<DG.DataFrame>
{
  return await _simulateBioreactor(initial, final, step,
    _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, 
    _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, 
    _yO2PInitial, _CYSTInitial, _VLInitial, _qinVal, _yO2inVal, 
    _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal);
}

//name: Bioreactor Demo
//description: In-browser simulation of complex phenomena.
//meta.demoPath: Bioreactors | Bioreactor
export async function demoScript(): Promise<any>  {
  const demoScript = new DemoScript('Bioreactor', 
    'No-code construction of complex phenomena simulators is provided by Datagrok WebAutosolver tool.'); 
  
  await demoScript    
    .step('Try', async () => 
    {
      const doeSimpleFunc: DG.Func = await grok.functions.eval('Bioreactors:BioreactorDemo');
      const doeSimpleFuncCall = doeSimpleFunc.prepare();
      
      const openModelFunc: DG.Func = await grok.functions.eval('Compute:openModelFromFuncall');
      const openModelFuncCall = openModelFunc.prepare({'funccall': doeSimpleFuncCall});
      openModelFuncCall.call();
    }, {description: 'Vary inputs and press "RUN".', delay: 0})
    //}, {description: 'Set inputs and press "RUN".', delay: 0})
    .step('Model', async () => {} , {description: 'Only declarative equations description is required.', delay: 0})    
    .step('Essence', async () => {} , {description: 'Simulation of controlled fab-arm exchange kinetic mechanism is performed here.', delay: 0})
    .step('Performance', async () => {} , {description: '1000 times faster than the previous version.', delay: 0})
    //.step('Complexity', async () => {} , {description: 'Each time you press "RUN", a system of 13 non-linear ordinary differential equations is solved.', delay: 0})
    .step('Complexity', async () => {} , {description: 'Each time you run computations, a system of 13 non-linear ordinary differential equations is solved.', delay: 0})
    .start();
}
