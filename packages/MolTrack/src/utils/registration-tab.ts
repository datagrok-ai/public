import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, Scope } from './constants';
import '../../css/moltrack.css';

let openedView: DG.ViewBase | null = null;

export class RegisrationView {
  view: DG.View;
  datasetInput: DG.InputBase | null = null;
  scopeInput: DG.InputBase | null = null;
  errorHandlingInput: DG.InputBase | null = null;
  mappingInput: DG.InputBase | null = null;
  gridDiv: HTMLDivElement;
  df: DG.DataFrame | null = null;

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register Entities';

    this.createInputs();
    this.gridDiv = ui.div('', 'moltrack-register-res-div');

    this.buildUI();
  }

  private createInputs() {
    this.datasetInput = ui.input.file('Dataset');
    this.scopeInput = this.createChoiceInput('Scope', Object.values(Scope));
    this.errorHandlingInput = this.createChoiceInput('Error handling', Object.values(ErrorHandling));
    this.mappingInput = ui.input.textArea('Mapping', { value: '' });
  }

  private buildUI() {
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerEntities());
    registerButton.classList.add('moltrack-run-register-button');

    const addToWorkspaceButton = ui.icons.add(() => {
      if (this.df)
        grok.shell.addTablePreview(this.df);
    }, 'Add registration results to workspace');

    this.view.setRibbonPanels([[addToWorkspaceButton]]);

    this.view.root.append(
      ui.divV([
        ui.divV([
          this.datasetInput?.root,
          this.scopeInput?.root,
          this.errorHandlingInput?.root,
          this.mappingInput?.root,
          registerButton,
        ], 'ui-form'),
        this.gridDiv,
      ]),
    );
  }

  private createChoiceInput<T>(label: string, choices: T[]): DG.ChoiceInput<T | null> {
    return ui.input.choice(label, {
      value: choices[0],
      items: choices,
    });
  }

  private async registerEntities() {
    ui.setUpdateIndicator(this.gridDiv, true);
    const file = this.datasetInput?.value;
    if (file && !file.id) {
      grok.shell.warning('Please upload a dataset.');
      ui.setUpdateIndicator(this.gridDiv, false);
      return;
    }

    try {
      this.df = await grok.functions.call('Moltrack:registerBulk', {
        csvFile: file,
        scope: this.scopeInput?.value,
        mapping: this.mappingInput?.value,
        errorHandling: this.errorHandlingInput?.value,
      });

      ui.empty(this.gridDiv);
      if (this.df) {
        this.df.name = `Register Entities: ${this.scopeInput?.value}`;
        this.gridDiv.append(this.df.plot.grid().root);
      }
    } catch (e: any) {
      grok.shell.error(`Failed to register entities: ${e.message}`);
    }

    ui.setUpdateIndicator(this.gridDiv, false);
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
