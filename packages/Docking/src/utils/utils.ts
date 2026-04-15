import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

import { BiostructureData } from '@datagrok-libraries/bio/src/pdb/types';

import { BINDING_ENERGY_COL, POSE_COL, setAffinity, setPose, TARGET_PATH } from './constants';

export function getFromPdb(pdbString: string): { [name: string]: number } {
  const result: { [name: string]: number } = {};
  const remarkRegex = /REMARK\s+\d+\s+([\w\s\(\)-]+)\.\s+([-\d.]+)/g;
  let match;
  while ((match = remarkRegex.exec(pdbString)) !== null)
    result[match[1].trim()] = parseFloat(match[2].trim());
  return result;
}

export function getFromPdbs(pdb: DG.SemanticValue): DG.DataFrame {
  const col = pdb.cell.column;
  const resultDf = DG.DataFrame.create(col.length);

  for (let idx = 0; idx < col.length; idx++) {
    const values = getFromPdb(col.get(idx));
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
  
  const receptorData: BiostructureData = {
    binary: false,
    data: receptor ? (await grok.dapi.files.readAsText(receptor)) : (await fetchPdbContent(receptorName.toUpperCase())),
    ext: receptor ? receptor.extension : 'pdb',
    options: {name: receptorName,},
  };
  return receptorData;
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
      const pdbContent = await response.text();
      return pdbContent;
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
  $(addColumnIcon)
    .css('color', '#2083d5')
    .css('position', 'absolute')
    .css('top', '2px')
    .css('left', '-12px')
    .css('margin-right', '5px');
  
  const idx = molecule.cell.rowIndex;
  const val: number = propertyCol.get(idx);
  const numStyle = {style: {width: '42px', textAlign: 'right', flexShrink: '0', flexGrow: '0'}};
  return ui.divH([addColumnIcon, ui.divText(val.toFixed(2), numStyle)],
    {style: {position: 'relative'}});
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
    const row = ui.div([
      ui.divText(cur.toFixed(2)),
      ui.divText('vs'),
      ui.divText(hov.toFixed(2)),
      ui.divText(`\u0394${sign}${delta.toFixed(2)}`),
    ]);
    row.style.cssText =
      'display:grid;grid-template-columns:42px 16px 42px auto;align-items:baseline;column-gap:12px;';
    const cells = row.children as HTMLCollectionOf<HTMLElement>;
    cells[0].style.cssText = 'text-align:right';
    cells[1].style.cssText = 'text-align:center;color:var(--grey-4)';
    cells[2].style.cssText = 'text-align:right';
    cells[3].style.cssText = `text-align:left;font-size:11px;${color ? 'color:' + color : ''}`;
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
  poseCol!.name = table.columns.getUnusedName(POSE_COL);
  setPose(poseCol!.name);
  const affinityCol = autodockResults.col(BINDING_ENERGY_COL);
  affinityCol!.name = table.columns.getUnusedName(BINDING_ENERGY_COL);
  setAffinity(affinityCol!.name);
  const processedTable = DG.DataFrame.fromColumns([poseCol!, affinityCol!]);
  affinityCol!.setTag(DG.TAGS.DESCRIPTION, affinityDescription);
  return processedTable;
}