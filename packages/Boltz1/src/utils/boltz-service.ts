import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as yaml from 'js-yaml';
import { BOLTZ_CONFIG_PATH, BOLTZ_PROPERTY_DESCRIPTIONS, BoltzResponse, Config } from './constants';
import { _package } from '../package';
import { getFromPdbs, prop } from './utils';

export class BoltzService {
  static async getBoltzConfigFolders(): Promise<string[]> {
    const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(BOLTZ_CONFIG_PATH, true);
    return targetsFiles
      .filter(folder => folder.isDirectory)
      .map(folder => folder.name);
  }

  static async runBoltz(config: string, msa: string): Promise<string> {
    const boltzContainer = await grok.dapi.docker.dockerContainers.filter('boltz').first();

    const body = {
      yaml: config,
      ...(msa && msa.trim() !== '' && { msa: msa }),
    };

    const params: RequestInit = {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body),
    };

    const response: Response = await grok.dapi.docker.dockerContainers.fetchProxy(boltzContainer.id, '/predict', params);
    const contentType = response.headers.get('content-type');
    const isJson = contentType === 'application/json';

    if (!response.ok) {
      if (isJson) {
        const json = await response.json();
        if (json['datagrok-error'])
          throw new Error(`Datagrok error: ${json['datagrok-error']}`);
        throw new Error(response.statusText);
      }
      const text = await response.text();
      throw new Error(`Error: ${text}`);
    }

    if (!isJson) {
      const text = await response.text();
      throw new Error(`Error: Boltz expected JSON response, got '${text}'.`);
    }

    const jsonResponse: BoltzResponse = await response.json();

    if (!jsonResponse.success) {
      _package.logger.error(jsonResponse.error);
      throw new Error('Prediction attempt failed: ' + jsonResponse.error);
    }

    return jsonResponse.result!;
  }

  static processBoltzResult(resultDf: DG.DataFrame) {
    const pdbCol = resultDf.columns.byName('pdb');
    const confidenceCol = resultDf.columns.byName('confidence_score');
    if (!pdbCol || !confidenceCol) return;

    pdbCol.semType = DG.SEMTYPE.MOLECULE3D;
    confidenceCol.meta.colors.setLinear([DG.Color.red, DG.Color.green]);
    confidenceCol.meta.format = '0.000';
    confidenceCol.setTag(DG.TAGS.DESCRIPTION, BOLTZ_PROPERTY_DESCRIPTIONS['confidence_score']);
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
        protein: { // Structure the sequence under the 'protein' field
          id: chainId, // Chain ID (e.g., 'A', 'B', etc.)
          sequence: sequence // The protein sequence
        }
      });

      const yamlString = yaml.dump(config);
      try {
        const result = DG.DataFrame.fromCsv(
          await grok.functions.call('Boltz1:runBoltz', { config: yamlString, msa: '' }));
        resultDf.append(result, true);
      } catch (err: any) {
        const msg = err?.message ?? String(err);
        _package.logger.error(`Folding failed for sequence ${index + 1}: ${msg}`);
        grok.shell.warning(`Skipped sequence ${index + 1}: ${msg}`);
      }
    }

    this.processBoltzResult(resultDf);
    await grok.data.detectSemanticTypes(resultDf);
    return resultDf;
  }  
  
  static async docking(df: DG.DataFrame, molecules: DG.Column, config: string): Promise<DG.DataFrame> {
    const configFile = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${config}`)).find((file) => file.extension === 'yaml')!;
    const msa = (await grok.dapi.files.list(`${BOLTZ_CONFIG_PATH}/${config}`)).find((file) => file.extension === 'a3m');

    let msaFile = '';
    if (msa)
      msaFile = await grok.dapi.files.readAsText(msa.fullPath);

    const existingConfig = yaml.load(await configFile.readAsString()) as any;
    let resultDf = DG.DataFrame.create();

    const isSmiles = molecules.meta.units === DG.UNITS.Molecule.SMILES;
    for (let [index, molecule] of molecules.toList().entries()) {
      const sequences = existingConfig.sequences;
      const constraints = existingConfig.constraints;

      const chainId = String.fromCharCode(65 + index);

      const ligandBlock = {
        ligand: {
          id: [chainId],
          smiles: isSmiles
            ? molecule
            : grok.chem.convert(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles),
        },
      };

      sequences.push(ligandBlock);
      constraints[0].pocket.binder = chainId;

      try {
        const updatedConfig = yaml.dump(existingConfig);
        const result = DG.DataFrame.fromCsv(
          await grok.functions.call('Boltz1:runBoltz', { config: updatedConfig, msa: msaFile }));
        resultDf.append(result, true);
      } catch (err: any) {
        const msg = err?.message ?? String(err);
        _package.logger.error(`Docking failed for ligand ${index + 1}: ${msg}`);
        grok.shell.warning(`Skipped ligand ${index + 1}: ${msg}`);
      } finally {
        sequences.pop();
      }
    }

    this.processBoltzResult(resultDf);
    await grok.data.detectSemanticTypes(resultDf);
    return resultDf;
  }
  
  static async boltzWidget(molecule: DG.SemanticValue): Promise<DG.Widget | null> {
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
      map[columnName] = prop(molecule, propertyCol!, result, BOLTZ_PROPERTY_DESCRIPTIONS);
    }
    result.appendChild(ui.tableFromMap(map));
    widget.root.append(result);
  
    return widget;
  }
}