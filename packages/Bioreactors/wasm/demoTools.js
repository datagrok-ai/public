// Custom demo tools

import {_simulateBioreactor} from './reactorAPI';

// custom demo
export async function customRun() {
    let initial = 0;
    let final = 1000;
    let step = 0.1;
    let _FFoxInitial = 0.268;
    let _KKoxInitial = 0.268;
    let _FFredInitial = 0;
    let _KKredInitial = 0;
    let _FfreeInitial = 0;
    let _KfreeInitial = 0;
    let _FKredInitial = 0;
    let _FKoxInitial = 0;
    let _MEAthiolInitial = 34;
    let _CO2Initial = 0.22;
    let _yO2PInitial = 0.209;
    let _CystamineInitial = 0; 
    let _VLInitial = 6.6;
    let _qinVal = 1;
    let _percentO2saturationVal = 100;
    let _yO2inVal = 0.209;
    let _pKa2MEAVal = 8.19;
    let _HVal = 1.072069378;
    let _TVal = 300;
    let _RVal = 0.082;
    let _PVal = 1;
    let _TimeToSwitchVal = 10;
  
    let outputNames = ['FFox(t)', 'KKox(t)', 'FFred(t)', 'KKred(t)', 
      'Ffree(t)', 'Kfree(t)', 'FKred(t)', 'FKox(t)', 'MEAthiol(t)', 
      'CO2(t)', 'yO2P(t)', 'Cystamine(t)', 'VL(t)'];
  
    let tablesCount = outputNames.length;
    let tables = Array(tablesCount);
      
    for (let i = 0; i < tablesCount; ++i) {
      tables[i] = DG.DataFrame.create(Math.trunc((final - initial) / step) + 1);
      tables[i].name = outputNames[i];
    }
  
    const count = 30;
  
    let view = grok.shell.newView('Simulation results');
  
    for (let i = 1; i <= count; ++i) {
      let output = await _simulateBioreactor(initial, final, step,
        _FFoxInitial, _KKoxInitial, _FFredInitial, _KKredInitial, _FfreeInitial, _KfreeInitial, 
        _FKredInitial, _FKoxInitial, _MEAthiolInitial, _CO2Initial, _yO2PInitial, _CystamineInitial, _VLInitial, 
        _qinVal, _percentO2saturationVal, _yO2inVal, _pKa2MEAVal, _HVal, _TVal, _RVal, _PVal,
        _TimeToSwitchVal + 5 * i);
      
      view.append(ui.h3(`Parameters set # ${i}`));
      view.append(ui.divH([
        DG.Viewer['lineChart'](output, { markerType: 'dot', sharex: 'true', multiAxis: 'true'}),
        DG.Viewer.grid(output)        
      ]));
  
      for (let j = 0; j < tablesCount; ++j) {
        if (!tables[j].columns.contains('t'))
          tables[j].columns.add(output.columns.byIndex(0));
  
        tables[j].columns.add(DG.Column.fromFloat32Array(
          `params #${i}`, output.columns.byIndex(j + 1).getRawData()));
      }
    }
    
    return tables;
  }
  
  // show result of custom run
  export function showCustomRunResults(tables) {
    let view = grok.shell.newView('Runs comparison');
  
    let count = tables.length;
  
    for (let i = count - 1; i >= 0; --i) {
      view.append(ui.h3('The function ' + tables[i].name));
      view.append(ui.divH([      
        DG.Viewer['lineChart'](tables[i], { markerType: 'dot', sharex: 'true', multiAxis: 'true'}),
        DG.Viewer.grid(tables[i])
      ]));
    }
  }
