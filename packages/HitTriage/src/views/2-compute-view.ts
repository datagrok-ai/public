import {HitTriageBaseView} from './hit-triage-base-view';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {HitTriageApp} from '../hit-triage-app';
import {chemDescriptorsDialog} from './descriptors-dialog';
import {ChemFunctionType, ChemPropNames} from './types';
import {ChemPropertyGroupMap, chemFunctionNameMap} from './consts';

/**
 * Enrichment of the molecular dataset.
 **/

const chemFunctionsMap: {[_ in ChemPropNames]: ChemFunctionType} = {
  'Descriptors': (t, c, p) => grok.chem.descriptors(t, c, p),
  'Toxicity Risks': getChemFunction('Toxicity Risks'),
  'Structural Alerts': getChemFunction('Structural Alerts'),
  'Chemical Properties': getChemFunction('Chemical Properties'),
};

function getChemFunction(funcName: Exclude<ChemPropNames, 'Descriptors'>): ChemFunctionType {
  const func = DG.Func.find({package: 'Chem', name: chemFunctionNameMap[funcName]})[0];
  return ((t: DG.DataFrame, col: string, props: string[]) => {
    const args = getChemFunctionArgs(t, col, props, funcName);
    return func.apply(args);
  });
}

function getChemFunctionArgs(
  t: DG.DataFrame, col: string, props: string[], funcName: ChemPropNames,
): {table: DG.DataFrame, molecules: string, [key: string]: any} {
  const args: {[_: string]: boolean} = {};
  const argNames = Object.values(ChemPropertyGroupMap).find((g) => g.name === funcName)!.values
    .map((v) => v.propertyName);
  argNames.forEach((argName) => args[argName] = false);
  props.forEach((p) => args[p] = true);
  return {
    table: t,
    molecules: col,
    ...args,
  };
}
export class ComputeView extends HitTriageBaseView {
  grid: DG.Grid;
  selectedDescriptors: {[key: string]: Array<string>} = {};

  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage | Compute';

    this.grid = this.template.hitsTable!.plot.grid();
    const content = ui.divV([
      ui.h1('This is where we enrich our data'),
      this.grid,
    ]);

    this.root.appendChild(content);
  }

  onActivated() {
    super.onActivated();
    this.process().then(() => {
    });
  }

  async process(): Promise<void> {
    return new Promise<void>(async (resolve) => {
      this.template.enrichedTable = this.template.hitsTable;

      chemDescriptorsDialog(async (descriptorsMap) => {
        Object.keys(descriptorsMap).forEach(async (key) => {
          await chemFunctionsMap[key as ChemPropNames](
            this.template.hitsTable!, this.template.hitsMolColumn, descriptorsMap[key]);
        });
        resolve();
      }, resolve);

      //const descriptors = ['HeavyAtomCount', 'NHOHCount'];
      //await grok.chem.descriptors(this.template.enrichedTable!, this.template.hitsMolColumn, descriptors);
      // another way to invoke calculations
      //await grok.functions.call('ChemDescriptors',
      // {table: session.enrichedDataFrame, smiles: session.sourceMoleculeColumn, descriptors: descriptors});

      this.grid.dataFrame = this.template.enrichedTable!;
    });
  }
}
