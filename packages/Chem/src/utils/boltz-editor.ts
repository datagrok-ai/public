import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BOLTZ_CONFIG_PATH, getBoltzConfigFolders } from '../package';

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
    const yaml = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${this.configInput.value}`)).find((file) => file.extension === 'yaml')!;
    const yamlContent = await yaml.readAsString();
    const match = yamlContent.match(/smiles:\s*([^\s]+)/);
    if (match) {
      const smilesInput = ui.input.molecule('Ligand', {value: match[1]});
      modelsSettingsDiv.appendChild(smilesInput.root);
    }
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