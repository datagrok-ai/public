import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {u2} from "@datagrok-libraries/utils/src/u2";
import * as yaml from 'js-yaml';

import { BoltzService } from '../utils/boltz-service';
import { _package } from '../package';

import '../css/boltz.css';

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
    this.divV.classList.add('boltz-app-wrapper');

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
    this.samplesInput.root.classList.add('boltz-samples-input');
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

  private createExampleButton(): HTMLElement {
    const exampleButton = ui.button('Example', () => this.fillProteinInputs());
    return exampleButton;
  }

  private fillProteinInputs(): void {
    this.setProteinInput(FUNC_PROPS_FIELD.PROTEIN + '_1', PROTEIN_1);
    this.setProteinInput(FUNC_PROPS_FIELD.PROTEIN + '_2', PROTEIN_2);
  }
  
  private setProteinInput(uniqueName: string, value: string): void {
    let input = this.sequenceInputs.get(uniqueName);
    if (!input) {
      const div = this.createInputSection(FUNC_PROPS_FIELD.PROTEIN, uniqueName);
      this.divV.insertBefore(div, this.samplesInput.root);
      input = this.sequenceInputs.get(uniqueName);
    }
    input!.value = value;
  }

  private getExampleYaml(): string {
    return yaml.dump({
      sequences: [
        { ENTITY_TYPE: FUNC_PROPS_FIELD.PROTEIN, id: FUNC_PROPS_FIELD.PROTEIN + '_1', sequence: PROTEIN_1, modifications: [] },
        { ENTITY_TYPE: FUNC_PROPS_FIELD.PROTEIN, id: FUNC_PROPS_FIELD.PROTEIN + '_2', sequence: PROTEIN_2, modifications: [] }
      ],
      constraints: []
    });
  }

  private createButtons(): void {
    const submitButton = ui.bigButton('Submit', async () => {
      const generatedYaml = this.generateYaml();
      const exampleYaml = this.getExampleYaml();
      if (generatedYaml === exampleYaml) {
        const cachedResults = DG.DataFrame.fromCsv(await _package.files.readAsText('demo_files/boltz_demo.csv'));
        grok.shell.addTableView(cachedResults);
        BoltzService.processBoltzResult(cachedResults);
        await grok.data.detectSemanticTypes(cachedResults);
      } else {
        // Run Boltz
      }
    });

    const addButton = this.createAddButton();
    const exampleButton = this.createExampleButton();

    Object.assign(submitButton.style, { margin: '0' });
    Object.assign(addButton.style, { margin: '0' });
    Object.assign(exampleButton.style, { margin: '0' });

    this.divV.append(
      ui.divH([addButton, exampleButton, submitButton], 'boltz-buttons-container')
    );
  }

  private getUniqueEntityName(type: FUNC_PROPS_FIELD): string {
    const count = (this.entityCounters.get(type) || 0) + 1;
    this.entityCounters.set(type, count);
    return `${type}_${count}`;
  }

  private createInputSection(type: FUNC_PROPS_FIELD, uniqueName: string): HTMLDivElement {
    const div = ui.divV([]);
    div.classList.add('boltz-input-section');

    const titleDiv = ui.divH([ui.divText(type)], 'boltz-title');

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
    tabControl.classList.add('boltz-tab-control');
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

    return ui.divH([binderChain.root, pocketChain.root, pocketResidues.root], { style: { justifyContent: 'space-between', width: '80%' } });
  }

  private createModifiableInput(uniqueName: string, div: HTMLDivElement): HTMLDivElement {
    const sequenceInp = ui.input.textArea('Sequence', {size: {width: 300, height: 50}});
    sequenceInp.root.classList.add('boltz-sequence-input');
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

  private getModificationsForSequence(sequenceName: string): { ccd: string; position: string }[] {
    return Array.from(this.modificationInputs.entries())
      .filter(([key]) => key.startsWith(sequenceName)) // Only modifications belonging to the sequence
      .map(([, inputs]) => {
        const position = inputs[0].value.trim();
        const ccd = inputs[1].value.trim();
        return position && ccd ? { position, ccd } : null;
      })
      .filter(Boolean) as { ccd: string; position: string }[];
  }

  private generateYaml(): string {
    const sequences = Array.from(this.sequenceInputs.entries())
      .filter(([name]) => !name.includes(" Position") && !name.includes(" CCD"))
      .map(([name, input]) => {
        const baseType = name.split('_')[0] as FUNC_PROPS_FIELD;
  
        if ([FUNC_PROPS_FIELD.PROTEIN, FUNC_PROPS_FIELD.RNA, FUNC_PROPS_FIELD.DNA].includes(baseType)) {
          return { 
            ENTITY_TYPE: baseType, 
            id: name, 
            sequence: input.value, 
            modifications: this.getModificationsForSequence(name)
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
    return yaml.dump({ sequences, constraints });
  }

  public getView(): DG.ViewBase {
    const view = DG.View.fromRoot(this.divV);
    view.name = 'Boltz-1';
    return view;
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

const PROTEIN_1 = "HHHHHHMAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV";
const PROTEIN_2 = "AVLFPVLLEVQRAVQLVVAPWLLQPAPGAPPPERLSLPAGRRLVVNLWATWCVLAALGLALGPPHRALVPGDPQAFRGLAHLVLGRHGNPLLLRVARAEGPRAVAAAFAGHPQLHLLPGSLALGQALELLLPWLQVLGPAPSCVLERGYLRRRQGAEGLPFRSWLREQGLAGLQARVARLLRALGPALAAELPPPVAAEALRAGLSAPVGGLLRAAARLLEAGQPPPGPEPLRAAALAVLAAIAAFLAAVAL"