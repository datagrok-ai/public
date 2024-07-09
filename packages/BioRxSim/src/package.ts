/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solve} from './solver';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: BioRxSim
//tags: model
//description: BioRx demo simulation
//input: dataframe inputs
//output: dataframe solution {caption: Bioreactor; viewer: Line chart(block: 100, segmentColumnName: "_Stage", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
//sidebar: @compute
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
export async function BioRxSim(inputs: DG.DataFrame): Promise<DG.DataFrame> {
  return await solve(inputs);
}
