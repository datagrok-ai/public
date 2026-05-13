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
).subscribe(async (molecule: string) => {
  if (!currentWidget)
    return;
  try {
    await updateRetrosynthesisWidget(molecule, currentWidget);
  } catch (e: any) {
    ui.empty(currentWidget.root);
    currentWidget.root.append(ui.divText(e?.message ?? 'An error occurred'));
  }
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

    const showWidget = (smiles: string) => {
      try {
        const widget = PackageFunctions.retroSynthesisPath(smiles);
        if (widget.root.parentElement !== retrosynthesisDiv) {
          ui.empty(retrosynthesisDiv);
          retrosynthesisDiv.append(widget.root);
        }
      } catch (e) {
        console.error(e);
        grok.shell.error('Invalid or empty molecule');
      }
    };

    showWidget(DEMO_MOLECULE);

    sketcher.onChanged.subscribe(() => {
      const smiles = sketcher.getSmiles();
      if (smiles)
        showWidget(smiles);
    });
    grok.shell.addPreview(view);
  }
}

