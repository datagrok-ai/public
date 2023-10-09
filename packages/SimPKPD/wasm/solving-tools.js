// Wasm-solving tools

import {callWasm} from './call-wasm';

export async function initSolvers() {
  await initOneCompartmentPK();
  await initTwoCompartmentPK();
}

export async function simulateOneCompartmentPK(initial, final, step,
  _depotInitial, _centrInitial, _periInitial, _effInitial, 
  _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val)
{
  const _tCount = Math.trunc((final - initial) / step) + 1;
  const _varsCount = 5;
  
  return callWasm(OneCompartmentPK, 'solveOneCompartmentPK',
    [ initial, final, step,
     _depotInitial, _centrInitial, _periInitial, _effInitial,
     _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val,
     _tCount, _varsCount ] );
}

export async function simulateTwoCompartmentPK(initial, final, step,
  _depotInitial, _centrInitial, _periInitial, _effInitial, 
  _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val)
{
  const _tCount = Math.trunc((final - initial) / step) + 1;
  const _varsCount = 5;
    
  return callWasm(TwoCompartmentPK, 'solveTwoCompartmentPK',
    [ initial, final, step,
     _depotInitial, _centrInitial, _periInitial, _effInitial,
     _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val,
     _tCount, _varsCount ] );
}
