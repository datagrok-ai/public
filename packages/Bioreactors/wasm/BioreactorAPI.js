import {callWasm} from '../wasm/callWasm';

export async function _initinitBioreactor() {
  await initBioreactor();
}

export async function _simulateBioreactor(initial, final, step,
  _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, 
  _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, 
  _yO2PInitial, _CYSTInitial, _VLInitial, _qinVal, _yO2inVal, 
  _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal)
{
  // these values are generated automatically by WebAutoSolver tools
  let _tCount = Math.trunc((final - initial) / step) + 1;
  let _varsCount = 14;
  
  return await callWasm(Bioreactor, 'solveBioreactor',
    [ initial, final, step,
     _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, 
     _KfreeInitial, _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, 
     _yO2PInitial, _CYSTInitial, _VLInitial,
     _qinVal, _yO2inVal, _HVal, _TVal, _RVal, _PVal, _TimeToSwitchVal,
     _tCount, _varsCount ] );
}
