/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {DEMO_MOLECULE} from './const';
import {updateRetrosynthesisWidget} from './utils';
import {Subject} from 'rxjs';
import {debounceTime} from 'rxjs/operators';
export * from './package.g';
export const _package = new DG.Package();

const moleculeStream$ = new Subject<string>();
let currentWidget: DG.Widget | null = null;

moleculeStream$.pipe(
  debounceTime(2000),
).subscribe({
  next: (molecule: string) => {
    if (currentWidget)
      updateRetrosynthesisWidget(molecule, currentWidget);
  },
  error: (error: any) => {
    if (currentWidget) {
      ui.empty(currentWidget.root);
      currentWidget.root.append(ui.divText(error?.message ?? 'An error occurred'));
    }
  },
});
export class PackageFunctions {
  @grok.decorators.panel({
    'meta': {'allowAddAsColumn': 'false', 'role': 'widgets', 'domain': 'chem'},
    'name': 'Chemistry | Retrosynthesis',
    'condition': 'true',
  })
  static retroSynthesisPath(
    @grok.decorators.param({'name': 'smiles', 'options': {'semType': 'Molecule'}}) molecule: string): DG.Widget {
    if (!currentWidget)
      currentWidget = new DG.Widget(ui.div('', 'retrosynthesis-widget-div'));
    ui.setUpdateIndicator(currentWidget.root, true, 'Calculating paths...');
    // Emit the new molecule value to the stream
    moleculeStream$.next(molecule);
    return currentWidget;
  }

  @grok.decorators.func()
  static retrosynthesisTopMenu(): void {
    (grok.shell.v as DG.TableView).addViewer('Retrosynthesis Viewer');
  }

  @grok.decorators.func({
    'meta': {
      'demoPath': 'Cheminformatics | Retrosynthesis',
    },
    'name': 'Retrosynthesis Demo',
    'description': 'Generate retrosynthesis paths',
  })
  static async retrosynthesisDemo(): Promise<void> {
    await grok.functions.call('Chem:initChemAutostart');
    const view = DG.View.create();
    view.name = 'Retrosynthesis Demo';

    const sketcher: DG.chem.Sketcher = new DG.chem.Sketcher();
    sketcher.setSmiles('COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5');
    const retrosynthesisDiv = ui.div('', 'retrosynthesis-demo');

    const container = ui.divH([
      sketcher.root,
      retrosynthesisDiv,
    ], {style: {height: '100%'}});

    view.append(container);

    let demoInited = false;
    sketcher.onChanged.subscribe(async () => {
      const smiles = sketcher.getSmiles();
      if (smiles) {
        try {
          ui.empty(retrosynthesisDiv);
          ui.setUpdateIndicator(retrosynthesisDiv, true, 'Calculating retrosyntehsis paths...');
          const widget = await PackageFunctions.retroSynthesisPath(!demoInited ? DEMO_MOLECULE : smiles);
          demoInited = true;
          retrosynthesisDiv.append(widget.root);
          ui.setUpdateIndicator(retrosynthesisDiv, false);
        } catch (e) {
          grok.shell.error('Invalid or empty molecule');
        }
      }
    });
    grok.shell.addPreview(view);
  }
}

