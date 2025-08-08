import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Scope } from '../utils/constants';
import dayjs from 'dayjs';

let openedView: DG.ViewBase | null = null;

export class RegistrationSingleView {
  view: DG.View;
  mainDiv: HTMLDivElement;
  compoundDiv: HTMLDivElement;
  batchDiv: HTMLDivElement;
  assayRunDiv: HTMLDivElement;
  assayResultsDiv: HTMLDivElement;
  gridDiv: HTMLDivElement;
  scopeInput: DG.InputBase | null = null;
  structureInput: DG.InputBase | null = null;
  epaBatchIdInput: DG.InputBase | null = null;
  assayNameInput: DG.InputBase | null = null;
  assayRunDateInput: DG.DateInput | null = null;
  inputTest;
  inputTest2;

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register Entity';
    this.createInputs();
    this.mainDiv = ui.div('', 'moltrack-register-single-div');
    this.compoundDiv = ui.div('', 'moltrack-register-single-div');
    this.batchDiv = ui.div('', 'moltrack-register-single-div');
    this.assayRunDiv = ui.div('', 'moltrack-register-single-div');
    this.assayResultsDiv = ui.div('', 'moltrack-register-single-div');
    this.gridDiv = ui.div('', 'moltrack-register-res-div');
    this.mainDiv.append(this.getScopeDiv(this.scopeInput!.value!));
    this.inputTest = ui.typeAhead('Country', {source: {
      local: ['USA', 'Ukraine', 'Antigua', 'United Kingdom', 'United Arab Emirates']}});
    const prop = DG.Property.fromOptions({
      'name': 'what',
      'inputType': 'Float',
      'min': 0,
      'max': 10,
      // @ts-ignore
      'showSlider': true,
      'step': 1,
    });
    this.inputTest2 = ui.input.forProperty(prop);
   this.buildUI();
  }

  private createInputs() {
    const scopeChoices = Object.values(Scope);
    this.scopeInput = ui.input.choice('Scope', {
      value: scopeChoices[0], items: scopeChoices, onValueChanged: () => {
        ui.empty(this.mainDiv);
        this.mainDiv.append(this.getScopeDiv(this.scopeInput!.value!));
      },
    });
    // Compound registration
    this.structureInput = ui.input.molecule('Structure');
    this.structureInput.classList.add('moltrack-register-single-input');

    // Batch registration
    this.epaBatchIdInput = ui.input.string('EPA Batch ID');

    // Assay run registration
    this.assayNameInput = ui.input.string('Assay run name');
    this.assayRunDateInput = ui.input.date('Assay run date', { value: dayjs() });
  }

  private getScopeDiv(scope: string): HTMLDivElement {
    let retVal = this.mainDiv;
    switch (scope) {
    case Scope.COMPOUNDS:
      ui.empty(this.compoundDiv);
      this.compoundDiv.append(this.scopeInput!.root);
      this.compoundDiv.append(this.structureInput!.root);
      retVal = this.compoundDiv;
      break;
    case Scope.BATCHES:
      ui.empty(this.batchDiv);
      this.batchDiv.append(this.scopeInput!.root);
      this.batchDiv.append(this.epaBatchIdInput!.root);
      this.batchDiv.append(this.structureInput!.root);
      retVal = this.batchDiv;
      break;
    case Scope.ASSAY_RUNS:
      ui.empty(this.assayRunDiv);
      this.assayRunDiv.append(this.scopeInput!.root);
      this.assayRunDiv.append(this.assayNameInput!.root);
      this.assayRunDiv.append(this.assayRunDateInput!.root);
      retVal = this.assayRunDiv;
      break;
    case Scope.ASSAY_RESULTS:
      ui.empty(this.assayRunDiv);
      this.assayResultsDiv.append(this.scopeInput!.root);
      this.assayResultsDiv.append(this.assayNameInput!.root);
      this.assayResultsDiv.append(this.assayRunDateInput!.root);
      this.assayResultsDiv.append(this.epaBatchIdInput!.root);
      this.assayResultsDiv.append(this.inputTest!.root);
      this.assayResultsDiv.append(this.inputTest2!.root);

      retVal = this.assayResultsDiv;
      break;
    }
    return retVal;
  }

  private buildUI() {
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerEntity());
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

  private async registerEntity() {
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
