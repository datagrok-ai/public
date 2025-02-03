import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BOLTZ_CONFIG_PATH, getBoltzConfigFolders } from '../package';
import $ from 'cash-dom';

export class BoltzBaseEditor {
  configInput!: DG.ChoiceInput<string>;
  modelsSettingsDiv = ui.inputs([]);

  constructor() {
    this.configInput = ui.input.choice('Configuration', {
      onValueChanged: async () =>  await this.createLigandsDiv(this.modelsSettingsDiv),
      nullable: false
    }) as DG.ChoiceInput<string>;

    this.initConfigs();
  }
  
  async createLigandsDiv(modelsSettingsDiv: HTMLElement): Promise<HTMLElement> {
    ui.empty(modelsSettingsDiv);
    const yaml = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${this.configInput.value}`))
        .find((file) => file.extension === 'yaml')!;
    const yamlContent = await yaml.readAsString();
    const match = yamlContent.match(/smiles:\s*([^\s]+)/);
    
    const ligandsContainer = ui.div([]);
  
    const addLigandInput = (smiles: string = '') => {
        const ligandDiv = ui.divH([]);
        ligandDiv.style.gap = '20px';
        ligandDiv.style.alignItems = 'center';

        const addButton = ui.iconFA('plus', () => addLigandInput(), 'Add ligand');
        const removeButton = ui.iconFA('minus', () => {
            ligandDiv.remove();
        }, 'Remove ligand');

        $(addButton).add(removeButton).css({
          'color': '#2083d5',
        });

        const smilesInput = ui.input.molecule('Ligand', { value: smiles });

        ligandDiv.append(smilesInput.root, addButton, removeButton);
        ligandsContainer.appendChild(ligandDiv);
    };

    if (match) {
        addLigandInput(match[1]);
    } else {
        addLigandInput();
    }

    modelsSettingsDiv.appendChild(ligandsContainer);
    return modelsSettingsDiv;
  }
  
  private async initConfigs(template?: string) {
    const templates = await getBoltzConfigFolders();
    this.configInput.items = templates;
    if (template)
      this.configInput.value = template;
  }
  
  public getEditor(): HTMLElement {
    return ui.div([
      this.configInput.root,
      this.modelsSettingsDiv,
    ], { style: { minWidth: '450px' }, classes: 'ui-form' });
  }
  
  public getParams() {
    return {
      config: this.configInput.value!
    };
  }
}