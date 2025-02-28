import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { _package } from '../package';
import {u2} from "@datagrok-libraries/utils/src/u2";

import * as yaml from 'js-yaml';

import '../css/docking.css';

export class Boltz1AppView {
  private sequenceInputs: Map<string, DG.InputBase> = new Map();
  private modificationInputs: Map<string, DG.InputBase[]> = new Map();
  private pocketInputs: Map<string, DG.InputBase> = new Map();
  private entityCounters: Map<FUNC_PROPS_FIELD, number> = new Map();
  private divV: HTMLDivElement;
  private menu!: DG.Menu;
  samplesInput!: DG.InputBase<number | null>;

  constructor() {
    this.divV = ui.divV([], 'ui-form');
    this.divV.classList.add('docking-app-wrapper');

    this.createHeader();
    this.createSamplesInput();
    this.createButtons();
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
    this.samplesInput = ui.input.int('Results', { value: 1 });
    this.samplesInput.root.classList.add('docking-samples-input');
    //this.inputs.set(this.samplesInput.caption, this.samplesInput);
    this.divV.append(this.samplesInput.root);
  }

  private createAddButton(): HTMLElement {
    this.menu = DG.Menu.popup();
    FUNC_PROPS_FIELDS.forEach(prop =>
      this.menu.item(prop, () => {
        const uniqueName = this.getUniqueEntityName(prop);
        const div = this.createInputSection(prop, uniqueName);
        this.divV.insertBefore(div, this.samplesInput.root);
        this.divV.scrollTop = this.divV.scrollHeight;
      })
    );

    const addButton = ui.button('Add', () => this.menu.show());
    return addButton;
  }

  private createButtons(): void {
    const submitButton = ui.bigButton('Submit', async () => {
      grok.shell.info('Generated YAML is logged to the console');
      console.log(this.generateYaml());
    });

    const addButton = this.createAddButton();

    Object.assign(submitButton.style, { margin: '0' });
    this.createAddButton();
    Object.assign(addButton.style, { margin: '0' });

    this.divV.append(
      ui.divH([addButton, submitButton], 'docking-buttons-container')
    );
  }

  private getUniqueEntityName(type: FUNC_PROPS_FIELD): string {
    const count = (this.entityCounters.get(type) || 0) + 1;
    this.entityCounters.set(type, count);
    return `${type}_${count}`;
  }

  private createInputSection(type: FUNC_PROPS_FIELD, uniqueName: string): HTMLDivElement {
    const div = ui.divV([]);
    div.classList.add('docking-input-section');

    const titleDiv = ui.divH([ui.divText(type)], 'docking-title');

    const deleteButton = ui.button(ui.icons.delete(() => {
      div.remove();
      this.sequenceInputs.delete(uniqueName);
      this.pocketInputs.delete(uniqueName);
    }), () => {});

    deleteButton.style.marginLeft = 'auto';
    titleDiv.append(deleteButton);
    div.append(titleDiv);
    
    switch (type) {
      case FUNC_PROPS_FIELD.LIGAND:
        div.append(this.createLigandTabControl(uniqueName));
        break;
      case FUNC_PROPS_FIELD.POCKET:
        div.append(this.createPocketRestraintInputs(uniqueName));
        break;
      default:
        div.append(this.createModifiableInput(uniqueName, div));
    }

    return div;
  }

  private createLigandTabControl(uniqueName: string): HTMLDivElement {
    const tabControl = ui.tabControl({
      'CCD': () => {
        const ccdInput = ui.input.string('CCD');
        this.sequenceInputs.set(`${uniqueName} ${ccdInput.caption}`, ccdInput);
        return ui.panel([ui.divH([ccdInput.root, ui.link('', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/', 'Find your ligand CCD')])]);
      },
      'Structure': () => {
        const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.INPLACE);
        const panel = ui.panel([sketcher.root]);
        panel.style.cssText = 'overflow: hidden !important; flex-grow: initial !important;';
        return panel;
      },
    }).root;
    tabControl.classList.add('docking-tab-control');
    return tabControl;
  }

  private createPocketRestraintInputs(uniqueName: string): HTMLElement {
    const binderChain = ui.input.choice('Binder Chain', { items: ['A', 'B', 'C', 'D'], nullable: true });
    const pocketChain = ui.input.choice('Pocket Chain', { items: ['A', 'B', 'C', 'D'], nullable: true });
    const pocketResidues = ui.input.list('Pocket Residues');

    pocketResidues.root.style.width = '150px';

    [binderChain, pocketChain].forEach(inp => inp.root.style.width = '50px');
    this.pocketInputs.set(`${uniqueName} ${binderChain.caption}`, binderChain);
    this.pocketInputs.set(`${uniqueName} ${pocketChain.caption}`, pocketChain);
    this.pocketInputs.set(`${uniqueName} ${pocketResidues.caption}`, pocketResidues);

    return ui.divH([binderChain.root, pocketChain.root, pocketResidues.root], { style: { justifyContent: 'stretch' } });
  }

  private createModifiableInput(uniqueName: string, div: HTMLDivElement): HTMLDivElement {
    const sequenceInp = ui.input.textArea('Sequence');
    sequenceInp.root.classList.add('docking-sequence-input');
    this.sequenceInputs.set(uniqueName, sequenceInp);

    let modificationCounter = 0;

    const addModificationsButton = ui.button('Add modifications', () => {
      ++modificationCounter;
      const position = ui.input.string('Position');
      const ccd = ui.input.string('CCD');
      ccd.root.append(ui.link('', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/', 'Find your ligand CCD'));

      div.append(ui.divH([position.root, ccd.root]));
      this.modificationInputs.set(`${uniqueName} ${modificationCounter}`, [position, ccd]);
    });

    return ui.divV([sequenceInp.root, addModificationsButton]);
  }

  private generateYaml(): string {
    const sequences = Array.from(this.sequenceInputs.entries())
      .filter(([name]) => !name.includes(" Position") && !name.includes(" CCD"))
      .map(([name, input]) => {
        const baseType = name.split('_')[0] as FUNC_PROPS_FIELD;

        const modifications = 
  
        if ([FUNC_PROPS_FIELD.PROTEIN, FUNC_PROPS_FIELD.RNA, FUNC_PROPS_FIELD.DNA].includes(baseType)) {
          return { 
            ENTITY_TYPE: baseType, 
            id: name, 
            sequence: input.value, 
            modifications: []
          };
        } else if (baseType === FUNC_PROPS_FIELD.LIGAND) {
          return { 
            ENTITY_TYPE: baseType, 
            id: name, 
            [input.value === 'SMILES' ? 'smiles' : 'ccd']: input.value 
          };
        }
      })
      .filter(Boolean);

    const constraints: any[] = [];

    this.pocketInputs.forEach((input, key) => {
      const uniqueName = key.split(' ')[0];
        
      const binderChain = this.pocketInputs.get(`${uniqueName} Binder Chain`)?.value;
      const pocketChain = this.pocketInputs.get(`${uniqueName} Pocket Chain`)?.value;
      const pocketResidues = this.pocketInputs.get(`${uniqueName} Pocket Residues`)?.value;
    
      const contacts = pocketResidues ? pocketResidues.map((resIdx: string) => {
        return [pocketChain, resIdx];
      }) : [];

      if (binderChain || pocketChain || pocketResidues) {
        constraints.push({
          pocket: {
            binder: binderChain,
            contacts: contacts
          }
        });
      }
    });
    console.log(sequences);
    console.log(constraints);
    return yaml.dump({ sequences, constraints });
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
  POCKET = 'Pocket',
}

export const FUNC_PROPS_FIELDS = [...Object.values(FUNC_PROPS_FIELD)];