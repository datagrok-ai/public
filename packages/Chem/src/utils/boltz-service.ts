import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as yaml from 'js-yaml';
import $ from 'cash-dom';
import { BOLTZ_CONFIG_PATH } from '../package';

export class BoltzService {
  static async getBoltzConfigFolders(): Promise<string[]> {
    const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(BOLTZ_CONFIG_PATH, true);
    return targetsFiles
      .filter(folder => folder.isDirectory)
      .map(folder => folder.name);
  }

  static async runBoltz(config: string, msa: string): Promise<string> {
    const container = await grok.dapi.docker.dockerContainers.filter('boltz').first();
    const body = {
      yaml: config,
      ...(msa && { msa: msa }),
    };

    const url = 'http://127.0.0.1:8001/predict';
    const response = await fetch(url, {
      method: 'POST',
      body: JSON.stringify(body),
      headers: {'Content-Type': 'application/json'},
    });

    const json = await response.json();
    return json['result'];
  }

  static async processBoltzResult(df: DG.DataFrame) {
    const pdbCol = df.columns.byName('pdb');
    const confidenceCol = df.columns.byName('confidence_score');
      
    pdbCol.semType = DG.SEMTYPE.MOLECULE3D;
    confidenceCol.meta.colors.setLinear([DG.Color.green, DG.Color.red]);
    confidenceCol.meta.format = '0.000';
    confidenceCol.setTag(DG.TAGS.DESCRIPTION, PROPERTY_DESCRIPTIONS['confidence_score']);
  }
  
  static async folding(df: DG.DataFrame, sequences: DG.Column): Promise<DG.DataFrame> {
    let resultDf = DG.DataFrame.create();
    
    for (let [index, sequence] of sequences.toList().entries()) {
      const config: Config = {
        version: 1,
        sequences: []
      };
  
      const chainId = String.fromCharCode(65 + index);
  
      config.sequences.push({
        protein: {  // Structure the sequence under the 'protein' field
          id: chainId,  // Chain ID (e.g., 'A', 'B', etc.)
          sequence: sequence  // The protein sequence
        }
      });
  
      const yamlString = yaml.dump(config);
      const result = DG.DataFrame.fromCsv(await grok.functions.call('Chem:runBoltz', { config: yamlString }));
      resultDf.append(result, true);
    }
  
    await this.processBoltzResult(resultDf);
    df.columns.add(resultDf.columns.byName('pdb'));
    df.columns.add(resultDf.columns.byName('confidence_score'));
  
    return df;
  }  
  
  static async docking(df: DG.DataFrame, molecules: DG.Column, config: string): Promise<DG.DataFrame> {
    const configFile = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${config}`)).find((file) => file.extension === 'yaml')!;
    const msa = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${config}`)).find((file) => file.extension === 'a3m');
    
    let msaFile;
    if (msa)
      msaFile = await grok.dapi.files.readAsText(msa.fullPath);
    
    const existingConfig = yaml.load(await configFile.readAsString()) as any;
    let resultDf = DG.DataFrame.create();

    for (let [index, molecule] of molecules.toList().entries()) {
      const sequences = existingConfig.sequences;
      const constraints = existingConfig.constraints;

      const chainId = String.fromCharCode(65 + index);
        
      const ligandBlock = {
        ligand: {
          id: [chainId],
          smiles: molecule
        }
      };
  
      sequences.push(ligandBlock);
      constraints[0].pocket.binder = chainId;
      const updatedConfig = yaml.dump(existingConfig);
      
      const result = DG.DataFrame.fromCsv(await grok.functions.call('Chem:runBoltz', { config: updatedConfig, msa: msaFile}));
      resultDf.append(result, true);
    }

    await this.processBoltzResult(resultDf);
    df.columns.add(resultDf.columns.byName('pdb'));
    df.columns.add(resultDf.columns.byName('confidence_score'));
  
    return df;
  }
  
  static async boltzWidget(molecule: DG.SemanticValue): Promise<DG.Widget<any> | null> {
    const value = molecule.value;
    const boltzResults: DG.DataFrame = getFromPdbs(molecule);
    const widget = new DG.Widget(ui.div([]));
  
    const targetViewer = await molecule.cell.dataFrame.plot.fromType('Biostructure', {
      pdb: value,
      zoom: true,
    });
    targetViewer.root.classList.add('bsv-container-info-panel');
    widget.root.append(targetViewer.root);
  
    const result = ui.div();
    const map: { [_: string]: any } = {};
    for (let i = 0; i < boltzResults!.columns.length; ++i) {
      const columnName = boltzResults!.columns.names()[i];
      const propertyCol = boltzResults!.col(columnName);
      map[columnName] = prop(molecule, propertyCol!, result);
    }
    result.appendChild(ui.tableFromMap(map));
    widget.root.append(result);
  
    return widget;
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

type SequenceModification = {
  position: number;  // Index of the residue, starting from 1
  ccd: string;       // CCD code of the modified residue
};

type ProteinEntity = {
  id: string | string[];  // Chain ID or multiple Chain IDs for identical entities
  sequence: string;       // Sequence (only for protein, dna, rna)
  msa?: string;            // MSA Path (only for protein)
  modifications?: SequenceModification[];  // List of modifications
};

type LigandEntity = {
  id: string | string[];  // Chain ID or multiple Chain IDs for identical entities
  smiles: string;         // SMILES (only for ligands)
  ccd: string;            // CCD (only for ligands)
};

type SequenceEntity = {
  protein?: ProteinEntity;  // Protein entity with sequence and msa
  ligand?: LigandEntity;    // Ligand entity with smiles and ccd
};

type Config = {
  version: number;
  sequences: SequenceEntity[];  // Array of sequence entities, each with an entity type
  constraints?: Constraint[];   // Optional constraints (can be empty)
};

type Constraint = {
  bond?: {
    atom1: [string, number, string];  // [CHAIN_ID, RES_IDX, ATOM_NAME]
    atom2: [string, number, string];  // [CHAIN_ID, RES_IDX, ATOM_NAME]
  };
  pocket?: {
    binder: string;                    // CHAIN_ID for the binder
    contacts: [string, number][];      // List of [CHAIN_ID, RES_IDX] for contacts
  };
};