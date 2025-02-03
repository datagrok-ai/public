import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BOLTZ_CONFIG_PATH } from '../package';
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

    const ligandsContainer = ui.div([]);
    let ligands: any[] = [];

    const ligandExists = (smiles: string) => {
      return yamlContent.includes(`smiles: ${smiles}`);
    };

    const updateYAML = async () => {
      const ligandsSection = ligands.map(ligand => {
        if (!ligandExists(ligand.smiles))
          return `  - ligand:\n      smiles: ${ligand.smiles}`;
        return '';
      }).join('\n');
      const updatedYAML = yamlContent.replace(/(ligand:\s*[\s\S]*?)(constraints:)/, `$1${ligandsSection}\n$2`);
      await grok.dapi.files.writeAsText(yaml, updatedYAML);
    };

    const removeFromYAML = async (smiles: string) => {
      const updatedYAML = yamlContent.replace(
        new RegExp(`-\\s*ligand:\\s*\\n\\s*smiles:\\s*${smiles}[\\s\\S]*?^(?=\\s*-\\s*ligand:|$)`, 'gm'), 
        ''
      );
    
      await grok.dapi.files.writeAsText(yaml, updatedYAML);
    }    

    const addLigandInput = (smiles: string = '') => {
      const ligandDiv = ui.divH([]);
      ligandDiv.style.gap = '20px';
      ligandDiv.style.alignItems = 'center';

      const addButton = ui.iconFA('plus', () => addLigandInput(), 'Add ligand');
      const removeButton = ui.iconFA('minus', () => {
        ligandDiv.remove();
        ligands = ligands.filter(l => l.smiles !== smilesInput.value);
        removeFromYAML(smilesInput.value);
      }, 'Remove ligand');

      $(addButton).add(removeButton).css({
        'color': '#2083d5',
      });

      const smilesInput = ui.input.molecule('Ligand', {
        value: smiles,
        onValueChanged: (value) => {
          const ligand = ligands.find(l => l.smiles === smiles);
          if (ligand) {
            ligand.smiles = value;
          }
          updateYAML();
        }
      });

      ligandDiv.append(smilesInput.root, addButton, removeButton);
      ligandsContainer.appendChild(ligandDiv);

      ligands.push({ smiles });
      updateYAML();
    };

    const initialLigands = yamlContent.match(/- ligand:\s*smiles:\s*([^\s]+)/g);
    if (initialLigands) {
      initialLigands.forEach(ligand => {
      const smiles = ligand.match(/smiles:\s*([^\s]+)/)?.[1];
      if (smiles) {
        if (!ligands.some(existingLigand => existingLigand.smiles === smiles)) {
          addLigandInput(smiles);
        }
      }
    });
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

export function prop(molecule: DG.SemanticValue, propertyCol: DG.Column, host: HTMLElement) : HTMLElement {
  const addColumnIcon = ui.iconFA('plus', () => {
    const df = molecule.cell.dataFrame;
    propertyCol.name = df.columns.getUnusedName(propertyCol.name);
    propertyCol.setTag(DG.TAGS.DESCRIPTION, PROPERTY_DESCRIPTIONS[propertyCol.name]);
    df.columns.add(propertyCol);
  }, `Calculate ${propertyCol.name} for the whole table`);

  ui.tools.setHoverVisibility(host, [addColumnIcon]);
  $(addColumnIcon)
    .css('color', '#2083d5')
    .css('position', 'absolute')
    .css('top', '2px')
    .css('left', '-12px')
    .css('margin-right', '5px');
  
  const idx = molecule.cell.rowIndex;
  return ui.divH([addColumnIcon, propertyCol.get(idx)], {style: {'position': 'relative'}});
}

export function getFromPdbs(pdb: DG.SemanticValue): DG.DataFrame {
  const col = pdb.cell.column;
  const resultDf = DG.DataFrame.create(col.length);
  
  for (let idx = 0; idx < col.length; idx++) {
    const pdbValue = col.get(idx);
    const remarkRegex = /REMARK\s+\d+\s+([^\d]+?)\s+([-\d.]+)/g;
    let match;
  
    while ((match = remarkRegex.exec(pdbValue)) !== null) {
      const colName = match[1].trim();
      const value = parseFloat(match[2].trim());
        
      let resultCol = resultDf.columns.byName(colName);
      if (!resultCol) resultCol = resultDf.columns.addNewFloat(colName);
        
      resultDf.set(colName, idx, value);
    }
  }

  return resultDf;
}

export async function getBoltzConfigFolders(): Promise<string[]> {
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(BOLTZ_CONFIG_PATH, true);
  return targetsFiles
    .filter(folder => folder.isDirectory)
    .map(folder => folder.name);
}

export const PROPERTY_DESCRIPTIONS: { [colName: string]: string } = {
  'confidence_score': 'Overall prediction quality score',
  'ptm': 'Global fold similarity measure',
  'iptm': 'Accuracy of chain interactions',
  'ligand_iptm': 'Confidence in ligand binding',
  'protein_iptm': 'Confidence in protein interactions',
  'complex_plddt': 'Average per-residue confidence score',
  'complex_iplddt': 'Confidence in chain interfaces',
  'complex_pde': 'Uncertainty in chain positioning',
  'complex_ipde': 'Uncertainty in interface docking',
  'chains_ptm': 'Confidence per individual chain',
  'pair_chains_iptm': 'Interaction accuracy between chains'
};
