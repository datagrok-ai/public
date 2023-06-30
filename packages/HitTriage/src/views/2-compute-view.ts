import {HitTriageBaseView} from './hit-triage-base-view';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {HitTriageApp} from '../hit-triage-app';
import {chemDescriptorsDialog} from './descriptors-dialog';

/**
 * Enrichment of the molecular dataset.
 **/

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
    const functions: {packageName: string, name: string}[] = [
      {packageName: 'Chem', name: 'addChemPropertiesColumns'},
      {packageName: 'Chem', name: 'addChemRisksColumns'},
      {packageName: 'Chem', name: 'structuralAlertsTopMenu'},
    ];
    return new Promise<void>(async (resolve) => {
      this.template.enrichedTable = this.template.hitsTable;

      chemDescriptorsDialog(async (resultMap) => {
        // await Promise.all(
        //   Object.keys(descriptorsMap).filter((key) => descriptorsMap[key] && descriptorsMap[key].length > 0)
        //     .map((key) => {
        //       return chemFunctionsMap[key as ChemPropNames](
        //     this.template.hitsTable!, this.template.hitsMolColumn, descriptorsMap[key]);
        //     }));
        const promises: Promise<any>[] = [];
        if (resultMap.descriptors && resultMap.descriptors.length > 0) {
          promises.push(
            grok.chem.descriptors(this.template.hitsTable!, this.template.hitsMolColumn, resultMap.descriptors,
            ));
        }
        Object.keys(resultMap.externals).forEach((funcName) => {
          const props = resultMap.externals[funcName];
          props['table'] = this.template.hitsTable!;
          props['molecules'] = this.template.hitsMolColumn;
          if (props)
            promises.push(grok.functions.call(funcName, props));
        });

        resolve();
      }, resolve, this.template.enrichedTable!, functions);

      //const descriptors = ['HeavyAtomCount', 'NHOHCount'];
      //await grok.chem.descriptors(this.template.enrichedTable!, this.template.hitsMolColumn, descriptors);
      // another way to invoke calculations
      //await grok.functions.call('ChemDescriptors',
      // {table: session.enrichedDataFrame, smiles: session.sourceMoleculeColumn, descriptors: descriptors});

      this.grid.dataFrame = this.template.enrichedTable!;
    });
  }
}
