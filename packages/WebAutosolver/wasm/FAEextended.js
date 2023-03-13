

var solveFAEextended = {
  arguments: {
    initial: {
      type: 'num'
    },
    final: {
      type: 'num'
    },
    step: {
      type: 'num'
    },
    _FFoxInitial: {
      type: 'num'
    },
    _KKoxInitial: {
      type: 'num'
    },
    _FFredInitial: {
      type: 'num'
    },
    _KKredInitial: {
      type: 'num'
    },
    _FfreeInitial: {
      type: 'num'
    },
    _KfreeInitial: {
      type: 'num'
    },
    _FKredInitial: {
      type: 'num'
    },
    _FKoxInitial: {
      type: 'num'
    },
    _MEAthiolInitial: {
      type: 'num'
    },
    _CO2Initial: {
      type: 'num'
    },
    _yO2PInitial: {
      type: 'num'
    },
    _CystamineInitial: {
      type: 'num'
    },
    _VLInitial: {
      type: 'num'
    },
    _qinVal: {
      type: 'num'
    },
    _percentO2saturationVal: {
      type: 'num'
    },
    _yO2inVal: {
      type: 'num'
    },
    _pKa2MEAVal: {
      type: 'num'
    },
    _HVal: {
      type: 'num'
    },
    _TVal: {
      type: 'num'
    },
    _RVal: {
      type: 'num'
    },
    _PVal: {
      type: 'num'
    },
    _TimeToSwitchVal: {
      type: 'num'
    },
    _tCount: {
      type: 'num'
    },
    _varsCount: {
      type: 'num'
    },
    _solution: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: '_tCount',
        value: 'data'
      },
      numOfColumns: {
        ref: '_varsCount',
        value: 'data'
      },
      names: ['t', 'FFox(t)', 'KKox(t)', 'FFred(t)', 'KKred(t)', 'Ffree(t)', 'Kfree(t)', 'FKred(t)', 'FKox(t)', 'MEAthiol(t)', 'CO2(t)', 'yO2P(t)', 'Cystamine(t)', 'VL(t)']
    }
  },
  output: {
    type: 'tableFromColumns',
    source: '_solution'
  }
}; // solveFAEextended

var FAEextended = undefined;

async function initFAEextended() {
  if (FAEextended === undefined) {
    console.log("Wasm not Loaded, Loading");
    FAEextended = await exportFAEextended();
    FAEextended.solveFAEextended = solveFAEextended;
  } else {
    console.log("Wasm Loaded, Passing");
  }
}