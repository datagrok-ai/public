import * as DG from 'datagrok-api/dg';


export class RulesForm {
  homoDimerInput: DG.InputBase<string>;
  heteroDimerInput: DG.InputBase<string>;

  // molSketcher: grok.chem.Sketcher;
  // monomerTypeInput: DG.ChoiceInput<string | null>;
  // polymerTypeInput: DG.ChoiceInput<string | null>;
  // monomerSymbolInput: DG.InputBase<string>;
  // monomerNameInput: DG.InputBase<string>;
  // monomerIdInput: DG.InputBase<number | null>;
  // monomerNaturalAnalogInput: DG.InputBase<string | null>;
  // rgroupsGrid: ItemsGrid;
  // metaGrid: ItemsGrid;
  // colors: {
  //   line: string,
  //   background: string,
  //   text: string
  // };
  // colorsEditor: ColorsEditor;
  // saveButton: HTMLButtonElement;
  // rgroupsGridRoot: HTMLElement;
  // private _molChanged: boolean = false;
  // get molChanged() { return this._molChanged; }
  // private saveValidationResult?: string | null = null;
  // private triggerMolChange: boolean = true; // makes sure that change is not triggered by copying the molecule from grid
  inputsTabControl: DG.TabControl;

  constructor() {
    this.homoDimerInput = ui.input.string('Homo dimer', {value: '($2)',
      onValueChanged: () => {}, nullable: false});
    this.heteroDimerInput = ui.input.string('Homo dimer', {value: '(#2)',
      onValueChanged: () => {}, nullable: false});

    const dimerInputsDiv = ui.divV([
      this.homoDimerInput,
      this.heteroDimerInput,
    ]);

    this.inputsTabControl = ui.tabControl({
      'Dimers': dimerInputsDiv
    }, false);
  }

  getForm() {
    this.inputsTabControl.root.classList.add('rules-manager-form-tab-control');
    this.inputsTabControl.header.style.marginBottom = '10px';

    return ui.divV([]);
    // return ui.divV([
    //   this.inputsTabControl.root,
    // ], {classes: 'ui-form', style: {paddingLeft: '10px', overflow: 'scroll'}});
  }
}
