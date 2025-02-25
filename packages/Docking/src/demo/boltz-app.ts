import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { _package } from '../package';
import {u2} from "@datagrok-libraries/utils/src/u2";

import * as yaml from 'js-yaml';

export class Boltz1AppView {
  private inputs: Array<{ name: string, value: DG.InputBase }> = [];
  private divV: HTMLDivElement;
  private submitButton!: HTMLButtonElement;
  private addIcon!: HTMLButtonElement;
  private menu!: DG.Menu;

  constructor() {
    this.divV = ui.divV([], 'ui-form');
    this.divV.style.overflow = 'unset';

    this.createHeader();
    this.createSamplesInput();
    this.createSubmitButton();
  }

  private createHeader(): void {
    const appHeader = u2.appHeader({
      iconPath: `${_package.webRoot}/images/boltz1_pred_figure.png`,
      learnMoreUrl: 'https://github.com/jwohlwend/boltz/blob/main/README.md',
      description: `- Seamless integration with Boltz.\n` +
        `- Predicts 3D structures for proteins, RNA, DNA, and small molecules.\n` +
        `- Supports modified residues, covalent ligands, glycans, and condition generation on pocket residues.\n` +
        `- Matches or outperforms the performance of AlphaFold3.`,
    });

    this.divV.append(appHeader);
  }

  private createSamplesInput(): void {
    const samplesInput = ui.input.int('Number of samples', { value: 1 });
    this.inputs.push({ name: 'Number of samples', value: samplesInput });
    this.createAddButton();
    this.divV.append(ui.divH([samplesInput.root, this.addIcon]));
  }

  private createAddButton(): void {
    this.menu = DG.Menu.popup();
    FUNC_PROPS_FIELDS.forEach((prop) => {
      this.menu.item(prop, () => {
        this.submitButton.remove();
        const div = this.createInputSection(prop);
        this.divV.append(div);
        this.divV.append(this.submitButton);
      });
    });

    this.addIcon = ui.button(ui.icons.add(() => this.menu.show()), () => {});
    this.addIcon.style.margin = '0px';
  }

  private createSubmitButton(): void {
    this.submitButton = ui.button('Submit', async () => {
      grok.shell.info('Generated YAML is logged to the console');
      this.inputs.forEach((input) => {
        console.log(this.generateYaml());
      });
    });
    this.submitButton.style.margin = '0px';

    this.divV.append(this.submitButton);
  }

  private createInputSection(name: string): HTMLDivElement {
    const div = ui.divV([]);
    div.style.border = '1px solid #ccc';
    div.style.padding = '10px';
    div.style.marginBottom = '10px';
    div.style.borderRadius = '8px';

    const titleDiv = ui.divH([]);
    const title = ui.divText(name);
    title.style.fontWeight = 'bold';
    title.style.fontSize = '14px';
    title.style.marginBottom = '10px';

    titleDiv.append(title);

    const deleteButton = ui.button(ui.icons.delete(() => {
      div.remove();
      this.inputs = this.inputs.filter((input) => input.name !== name);
    }), () => {});
    deleteButton.style.marginLeft = 'auto';
    titleDiv.append(deleteButton);
    div.append(titleDiv);

    const sequenceInp = ui.input.string(name);
    this.inputs.push({ name, value: sequenceInp });

    if (name === 'Ligand') {
      div.append(this.createLigandTabControl());
    } else if (name === 'Pocket Restraint') {
      div.append(...this.createPocketRestraintInputs());
    } else {
      div.append(this.createModifiableInput(sequenceInp, div));
    }

    return div;
  }

  private createLigandTabControl(): HTMLDivElement {
    const tabControl =  ui.tabControl({
      'CCD': () => {
        const ccdInput = ui.input.string('CCD');
        this.inputs.push({ name: 'Ligand CCD', value: ccdInput });
        return ui.panel([ui.divH([ccdInput.root, ui.link('', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/', 'Find your ligand CCD')])]);
      },
      'SMILES': () => {
        const smilesInput = ui.input.molecule('SMILES');
        this.inputs.push({ name: 'Ligand SMILES', value: smilesInput });
        return ui.panel([smilesInput.root]);
      },
    }).root;
    tabControl.style.height = '150px';
    return tabControl;
  }

  private createPocketRestraintInputs(): HTMLElement[] {
    const binderChainChoice = ui.input.choice('Binder Chain', { items: ['A', 'B', 'C', 'D'] });
    const pocketChainChoice = ui.input.choice('Pocket Chain', { items: ['A', 'B', 'C', 'D'] });
    const pocketResiduesList = ui.input.string('Pocket Residues (comma separated)');

    this.inputs.push({ name: 'Binder Chain', value: binderChainChoice });
    this.inputs.push({ name: 'Pocket Chain', value: pocketChainChoice });
    this.inputs.push({ name: 'Pocket Residues', value: pocketResiduesList });

    return [binderChainChoice.root, pocketChainChoice.root, pocketResiduesList.root];
  }

  private generateYaml(): string {
    const sequences = this.inputs
      .map(input => {
        const { name: inputName, value: { value } } = input;
        let sequenceObj: any;
  
        if ([FUNC_PROPS_FIELD.PROTEIN, FUNC_PROPS_FIELD.RNA, FUNC_PROPS_FIELD.DNA].includes(inputName as FUNC_PROPS_FIELD)) {
          const modifications = this.inputs.filter(i => i.name.includes(value));
          sequenceObj = {
            ENTITY_TYPE: inputName,
            id: value,
            sequence: value,
            position: modifications.find(m => m.name.includes('Position'))?.value.value,
            ccd: modifications.find(m => m.name.includes('CCD'))?.value.value,
          };
        }
        else if (inputName.includes(FUNC_PROPS_FIELD.LIGAND)) {
          const ligandType = value === 'SMILES' ? 'smiles' : 'ccd';
          sequenceObj = {
            ENTITY_TYPE: inputName,
            id: value,
            [ligandType]: value,
          };
        }
  
        return sequenceObj;
      })
      .filter(Boolean);
  
    const yamlData = {
      sequences,
      constraints: [],
    };
  
    const yamlString = yaml.dump(yamlData);
    return yamlString;
  }  

  private createModifiableInput(sequenceInp: DG.InputBase, div: HTMLDivElement): HTMLDivElement {
    const modificationsDiv = ui.divV([]);
    modificationsDiv.style.display = 'none';
    let modificationsAdded = false;

    const addModificationsButton = ui.button(ui.icons.settings(() => {
      if (!modificationsAdded) {
        const position = ui.input.string('Position');
        const ccd = ui.input.string('CCD');
        ccd.root.append(ui.link('', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/', 'Find your ligand CCD'));
        modificationsDiv.append(position.root, ccd.root);
        this.inputs.push({ name: `${sequenceInp.stringValue} Position`, value: position });
        this.inputs.push({ name: `${sequenceInp.stringValue} CCD`, value: ccd });
        div.append(modificationsDiv);
        modificationsAdded = true;
      }

      modificationsDiv.style.display = modificationsDiv.style.display === 'none' ? 'block' : 'none';
    }), () => {}, 'Add modifications');

    addModificationsButton.style.margin = '0px';

    return ui.divH([sequenceInp.root, addModificationsButton]);
  }

  public getView(): DG.ViewBase {
    return DG.View.fromRoot(this.divV);
  }
}

export enum FUNC_PROPS_FIELD {
  PROTEIN = 'Protein',
  RNA = 'RNA',
  DNA = 'DNA',
  LIGAND = 'Ligand',
  POCKET_RESTRAINT = 'Pocket Restraint',
}

export const FUNC_PROPS_FIELDS = [...Object.values(FUNC_PROPS_FIELD)];