// sensitivityAnalysis.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, getVariedNumericalInputColumns} from './inputTools';
import {getOutputColumns} from './outputTools';
import {checkSize} from './utils';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

export class VarianceBasedSenstivityAnalysis {
  private samplesCount: number;
  private runsCount: number;
  private fixedInputs: FixedInputItem[];
 
  private variedInputs: VariedNumericalInputValues[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[];  
  
  constructor(func: DG.Func, fixedInputs: FixedInputItem[],
    variedInputs: VariedNumericalInputInfo[], samplesCount: number) {
    // check size
    checkSize(samplesCount);

    this.samplesCount = samplesCount;
    this.fixedInputs = fixedInputs;
    
    this.func = func;

    const numericalColumns = getVariedNumericalInputColumns(samplesCount, variedInputs);

    const dimension = variedInputs.length;   
    
    this.variedInputs = [...Array(dimension).keys()].map(i => 
      ({prop: variedInputs[i].prop, 
        min: variedInputs[i].min, 
        max: variedInputs[i].max,
        column: numericalColumns[i]})
    );

    this.runsCount = samplesCount * (dimension + 2);

    this.funcCalls = [];

    for (let i = 0; i < this.runsCount; ++i) {
      let inputs: any = {};
    
      for (const input of fixedInputs)
        inputs[input.name] = input.value;
    
      for (const input of this.variedInputs)
        inputs[input.prop.name] = input.column.get(i);
          
      this.funcCalls.push(func.prepare(inputs));
    }
  }

  private async run(): Promise<void> {    
    await Promise.all(this.funcCalls.map((call) => call.call()));
  }

  async perform(): Promise<DG.DataFrame> {
    await this.run();

    const columns = this.variedInputs.map((varInput) => varInput.column as DG.Column)
    
    const table = DG.DataFrame.fromColumns(columns);
    table.name = `Sensitivity Analysis of ${this.func.friendlyName}`;    

    for (const col of getOutputColumns(this.funcCalls))
      table.columns.add(col);

    return table;
  }
};