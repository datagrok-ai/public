import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let openedView: DG.ViewBase | null = null;

export class RegistrationCompoundView {
  view: DG.View;
  gridDiv: HTMLDivElement;
  molfileInput: DG.InputBase | null = null;
  smilesInput: DG.InputBase | null = null;
  sketcherWidget: HTMLElement | null = null;

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register compound';
    this.createInputs();
    this.gridDiv = ui.div('', 'moltrack-register-res-div');
    this.buildUI();
  }

  private createInputs() {
    this.sketcherWidget = grok.chem.sketcher(this.onChange.bind(this), 'CC(=O)Oc1ccccc1C(=O)O)');
    this.molfileInput = ui.input.textArea('',{value: ""});
    this.smilesInput = ui.input.string('',{value: ""});
  }

  private buildUI() {
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerCompound());
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
        ]),
          registerButton,
        ], 'ui-form'),
        this.gridDiv,
      ]),
    );
  }

  private async registerCompound() {
    console.log("SMILES: ",this.smilesInput!.value);
  }

  private onChange(smiles: string, molfile: string) {
    if (this.smilesInput){
      this.smilesInput!.value = smiles;
      this.molfileInput!.value = molfile;
    }
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
