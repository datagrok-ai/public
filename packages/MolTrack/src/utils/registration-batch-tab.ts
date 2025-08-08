import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { createModifiableInput } from "./utils";

let openedView: DG.ViewBase | null = null;

export class RegistrationBatchView {
  view: DG.View;
  gridDiv: HTMLDivElement;
  propertiesDiv: HTMLDivElement;
  molfileInput: DG.InputBase | null = null;
  smilesInput: DG.InputBase | null = null;
  sketcherWidget: HTMLElement | null = null;
  epaBatchIdInput: DG.InputBase | null = null;
  modificationInputs: Map<string, any> = new Map();
  properties = ["batch_detail_1", "batch_detail_2"];

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register batch';
    this.createInputs();
    this.gridDiv = ui.div('', 'moltrack-register-res-div');
    this.propertiesDiv = ui.div('', 'moltrack-register-res-div');
    this.buildUI();
  }

  private createInputs() {
    this.sketcherWidget = grok.chem.sketcher(this.onChange.bind(this), 'CC(=O)Oc1ccccc1C(=O)O)');
    this.molfileInput = ui.input.textArea('', { value: '' });
    this.smilesInput = ui.input.string('', { value: '' });
    this.epaBatchIdInput = ui.input.string('EPA Batch ID');
  }

  private buildUI() {
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerBatch());
    registerButton.classList.add('moltrack-run-register-button');

    const addToWorkspaceButton = ui.icons.add(() => {
      // if (this.df)
      //     grok.shell.addTablePreview(this.df);
    }, 'Add registration results to workspace');

    this.view.setRibbonPanels([[addToWorkspaceButton]]);
    this.view.root.append(
      ui.divV([
        ui.divV([
          ui.divV([
            this.sketcherWidget,
            this.epaBatchIdInput!.root,
            createModifiableInput(this.propertiesDiv, this.modificationInputs, this.properties),
            this.propertiesDiv,
          ]),

          registerButton,
        ], 'ui-form'),
        this.gridDiv,
      ]),
    );
  }

  private onChange(smiles: string, molfile: string) {
    if (this.smilesInput) {
      this.smilesInput!.value = smiles;
      this.molfileInput!.value = molfile;
    }
  }

  private async registerBatch() {
    console.log('SMILES: ', this.smilesInput!.value);
    const results: { property: string, value: string }[] = [];

    for (const [, [propertyInput, valueInput]] of this.modificationInputs) {
      const property = propertyInput.value?.trim();
      const value = valueInput.value?.trim();

      if (property && value) {
        results.push({ property, value });
      }
    }
    console.log("EPA: ",this.epaBatchIdInput!.value);
    console.log("Details: ", results);
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
