import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BiostructureData } from '@datagrok-libraries/bio/src/pdb/types';

import { BINDING_ENERGY_COL, POSE_COL, setAffinity, setPose, TARGET_PATH } from './constants';

export function getRemarksFromPdb(pdbString: string): { [name: string]: number } {
  const result: { [name: string]: number } = {};
  const remarkRegex = /REMARK\s+\d+\s+([\w\s\(\)-]+)\.\s+([-\d.]+)/g;
  let match;
  while ((match = remarkRegex.exec(pdbString)) !== null)
    result[match[1].trim()] = parseFloat(match[2].trim());
  return result;
}

export function getRemarksFromPdbs(pdb: DG.SemanticValue): DG.DataFrame {
  const col = pdb.cell.column;
  const resultDf = DG.DataFrame.create(col.length);

  for (let idx = 0; idx < col.length; idx++) {
    const values = getRemarksFromPdb(col.get(idx));
    for (const [colName, value] of Object.entries(values)) {
      let resultCol = resultDf.columns.byName(colName);
      if (!resultCol) resultCol = resultDf.columns.addNewFloat(colName);
      resultDf.set(colName, idx, value);
    }
  }

  return resultDf;
}
  
export async function getReceptorData(pdb: string): Promise<BiostructureData> {
  const match = pdb.match(/REMARK\s+1\s+receptor\.\s+(\S+)\s+J\./);
  const receptorName = match ? match[1] : '';
  const receptor = await getReceptorFile(receptorName);
  
  return {
    binary: false,
    data: receptor ? (await grok.dapi.files.readAsText(receptor)) : (await fetchPdbContent(receptorName.toUpperCase())),
    ext: receptor ? receptor.extension : 'pdb',
    options: {name: receptorName,},
  } as BiostructureData;
}
  
export async function getReceptorFile(receptorName: string): Promise<DG.FileInfo | undefined> {
  const receptorPath = `${TARGET_PATH}/${receptorName}`;
  const folderExists = await grok.dapi.files.exists(receptorPath);
  if (!folderExists) return undefined;
  
  const files = await grok.dapi.files.list(receptorPath);
  return files.find(file => ['pdbqt', 'pdb'].includes(file.extension));
}
  
async function fetchPdbContent(pdbId: string, format: string = 'pdb'): Promise<string> {
  const baseUrl = 'https://files.rcsb.org/download';
  const url = `${baseUrl}/${pdbId}.${format}`;
  
  try {
    const response = await grok.dapi.fetchProxy(url);
    if (response.ok) {
      return  await response.text();
    }
  } catch (error) {

  }
  return '';
}

export function prop(molecule: DG.SemanticValue, propertyCol: DG.Column, host: HTMLElement, descriptions: { [colName: string]: string }) : HTMLElement {
  const addColumnIcon = ui.iconFA('plus', () => {
    const df = molecule.cell.dataFrame;
    propertyCol.name = df.columns.getUnusedName(propertyCol.name);
    propertyCol.setTag(DG.TAGS.DESCRIPTION, descriptions[propertyCol.name]);
    df.columns.add(propertyCol);
  }, `Calculate ${propertyCol.name} for the whole table`);

  ui.tools.setHoverVisibility(host, [addColumnIcon]);
  addColumnIcon.classList.add('docking-add-property-icon');

  const idx = molecule.cell.rowIndex;
  const val: number = propertyCol.get(idx);
  const numEl = ui.divText(val.toFixed(2), 'docking-property-value');
  const wrapper = ui.divH([addColumnIcon, numEl], 'docking-property-cell');
  return wrapper;
}

export function buildComparisonTable(
  currentValues: { [name: string]: number },
  hoveredValues: { [name: string]: number },
): HTMLElement {
  const map: { [_: string]: any } = {};
  for (const name of Object.keys(currentValues)) {
    const cur = currentValues[name];
    const hov = hoveredValues[name];
    if (hov === undefined) {
      map[name] = ui.divText(cur.toFixed(2));
      continue;
    }
    const delta = hov - cur;
    const sign = delta >= 0 ? '+' : '';
    const color = delta < 0 ? 'var(--green-2)' : delta > 0 ? 'var(--red-2)' : '';
    const deltaEl = ui.divText(`\u0394${sign}${delta.toFixed(2)}`, 'docking-score-delta');
    if (color) deltaEl.style.color = color;
    const row = ui.div([
      ui.divText(cur.toFixed(2), 'docking-score-current'),
      ui.divText('vs', 'docking-score-separator'),
      ui.divText(hov.toFixed(2), 'docking-score-hovered'),
      deltaEl,
    ], 'docking-comparison-row');
    map[name] = row;
  }
  return ui.tableFromMap(map);
}

export function formatColumns(autodockResults: DG.DataFrame) {
  for (let col of autodockResults.columns.numerical)
    col.meta.format = '0.00';
}

export function addColorCoding(column: DG.GridColumn) {
  column.isTextColorCoded = true;
  column.column!.meta.colors.setLinear([DG.Color.green, DG.Color.red]);
}

export function processAutodockResults(autodockResults: DG.DataFrame, table: DG.DataFrame): DG.DataFrame {
  const affinityDescription = 'Estimated Free Energy of Binding.\
    Lower values correspond to stronger binding.';
  const poseCol = autodockResults.col(POSE_COL);
  const affinityCol = autodockResults.col(BINDING_ENERGY_COL);
  // Both columns are required. Missing either means AutoDock didn't
  // produce the expected output (e.g. the docking-autodock container
  // was stopped, the target folder was unreadable, or the pipeline
  // returned partial results). Surface this with a shell warning and
  // return an empty DataFrame instead of crashing on `null!.name`.
  if (!poseCol || !affinityCol) {
    // `auto-dock-app.ts#getAutodockResults` creates a fallback DataFrame
    // with only a `pose` STRING column (no `binding energy`) when every
    // ligand failed to dock. The per-ligand error messages are written
    // into that pose column. Surface them so the user sees the real
    // problem instead of a generic "missing column" notice.
    if (poseCol && !affinityCol && poseCol.type === DG.TYPE.STRING) {
      const errorMessages: string[] = [];
      for (let i = 0; i < poseCol.length && errorMessages.length < 3; i++) {
        const v = poseCol.get(i);
        if (typeof v === 'string' && v.trim().length > 0)
          errorMessages.push(v.trim().slice(0, 200));
      }
      const msgSuffix = errorMessages.length > 0
        ? ` Docking errors: ${errorMessages.join(' | ')}`
        : ' Docker container reachable but no valid poses produced — check container logs.';
      grok.shell.warning(`AutoDock failed to produce binding-energy results.${msgSuffix}`);
      return DG.DataFrame.create();
    }
    const missing = [
      !poseCol ? POSE_COL : null,
      !affinityCol ? BINDING_ENERGY_COL : null,
    ].filter(Boolean).join(', ');
    const actualCols = autodockResults.columns.names().join(', ');
    grok.shell.warning(
      `AutoDock results are missing expected column(s): ${missing}. ` +
      `Actual columns returned: [${actualCols}]. ` +
      `Check that the docking-autodock container is running and the target folder is accessible.`);
    return DG.DataFrame.create();
  }
  poseCol.name = table.columns.getUnusedName(POSE_COL);
  setPose(poseCol.name);
  affinityCol.name = table.columns.getUnusedName(BINDING_ENERGY_COL);
  setAffinity(affinityCol.name);
  const processedTable = DG.DataFrame.fromColumns([poseCol, affinityCol]);
  affinityCol.setTag(DG.TAGS.DESCRIPTION, affinityDescription);
  return processedTable;
}