// sensitivityAnalysis.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, createVariedNumericalInputColumns} from './inputTools';
import {getOutputColumns} from './outputTools';
import {checkSize} from './utils';

export class VarianceBasedSenstivityAnalysis {
  private samplesCount: number;
  private runsCount: number;
  private fixedInputs: FixedInputItem[];
  private variedInputs: VariedNumericalInputInfo[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[];  
  
  constructor(func: DG.Func, fixedInputs: FixedInputItem[],
    variedInputs: VariedNumericalInputInfo[], samplesCount: number) {
    // check size
    checkSize(samplesCount);

    this.samplesCount = samplesCount;
    this.fixedInputs = fixedInputs;
    this.variedInputs = variedInputs;

    this.func = func;
    
    createVariedNumericalInputColumns(samplesCount, this.variedInputs);

    this.runsCount = samplesCount * (variedInputs.length + 2);

    this.funcCalls = [];

    for (let i = 0; i < this.runsCount; ++i) {
      let inputs: any = {};
    
      for (const input of fixedInputs)
        inputs[input.name] = input.value;
    
      for (const input of variedInputs)
        inputs[input.name] = input.column?.get(i);
          
      this.funcCalls.push(func.prepare(inputs));
    }
  }

  private async run(): Promise<void> {
    for (const funcCall of this.funcCalls)
      await funcCall.call();
  }

  async perform(): Promise<DG.DataFrame> {
    await this.run();

    const columns: DG.Column[] = [];

    for (const input of this.variedInputs)
      columns.push(input.column as DG.Column);
    
    const table = DG.DataFrame.fromColumns(columns);
    table.name = `Sensitivity Analysis of ${this.func.friendlyName}`;    

    for (const col of getOutputColumns(this.funcCalls))
      table.columns.add(col);

    return table;
  }
};