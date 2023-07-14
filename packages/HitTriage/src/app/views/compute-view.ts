import {HitTriageBaseView} from './base-view';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {HitTriageApp} from '../hit-triage-app';
import {chemFunctionsDialog} from '../dialogs/functions-dialog';

/**
 * Enrichment of the molecular dataset.
 **/

export class ComputeView extends HitTriageBaseView {
  grid: DG.Grid;
  selectedDescriptors: {[key: string]: Array<string>} = {};

  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Compute';

    this.grid = this.app.dataFrame?.plot.grid()!;
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
    const functions: {packageName: string, name: string}[] =
      this.app.template!.compute.functions.map((func) =>({packageName: func.package, name: func.name}));
    return new Promise<void>(async (resolve) => {
      chemFunctionsDialog(async (resultMap) => {
        const promises: Promise<any>[] = [];
        if (resultMap.descriptors && resultMap.descriptors.length > 0) {
          promises.push(
            grok.chem.descriptors(this.app.dataFrame!, this.app.molColName!, resultMap.descriptors,
            ));
        }
        Object.keys(resultMap.externals).forEach((funcName) => {
          const props = resultMap.externals[funcName];
          props['table'] = this.app.dataFrame!;
          props['molecules'] = this.app.molColName!;
          if (props)
            promises.push(grok.functions.call(funcName, props));
        });
        await Promise.all(promises);
        resolve();
      }, resolve, functions, this.app.template?.compute.descriptors.enabled);
    });
  }
}
