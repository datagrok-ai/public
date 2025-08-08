import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

let openedView: DG.ViewBase | null = null;

export class RegistrationAssayResultView {
  view: DG.View;
  mainDiv: HTMLDivElement;
  gridDiv: HTMLDivElement;
  epaBatchIdInput: DG.InputBase | null = null;
  assayNameInput: DG.InputBase | null = null;
  assayRunDateInput: DG.DateInput | null = null;

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register assay result';
    this.createInputs();
    this.mainDiv = ui.div('', 'moltrack-register-single-div');
    this.gridDiv = ui.div('', 'moltrack-register-res-div');
    this.buildUI();
  }

  private createInputs() {
    this.epaBatchIdInput = ui.input.string('EPA Batch ID');
    this.assayNameInput = ui.input.string('Assay run name');
    this.assayRunDateInput = ui.input.date('Assay run date', { value: dayjs() });
  }

  private buildUI() {
    this.mainDiv.append(this.assayNameInput!.root);
    this.mainDiv.append(this.assayRunDateInput!.root);
    this.mainDiv.append(this.epaBatchIdInput!.root);
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerAsayResult());
    registerButton.classList.add('moltrack-run-register-button');

    const addToWorkspaceButton = ui.icons.add(() => {
      // if (this.df)
      //     grok.shell.addTablePreview(this.df);
    }, 'Add registration results to workspace');

    this.view.setRibbonPanels([[addToWorkspaceButton]]);
    this.view.root.append(
      ui.divV([
        ui.divV([
          this.mainDiv,
          registerButton,
        ], 'ui-form'),
        this.gridDiv,
      ]),
    );
  }

  private async registerAsayResult() {
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
