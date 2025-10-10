import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { prop } from '@datagrok-libraries/utils/src/add-icon-utils';
import { CACHED_RESULTS } from './constants';

export function getFromPdbs(pdb: DG.SemanticValue): DG.DataFrame {
  const col = pdb.cell.column;
  if (CACHED_RESULTS.has(col.toString()))
    return CACHED_RESULTS.get(col.toString())!;
      
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
    
  CACHED_RESULTS.set(col.toString(), resultDf);
    
  return resultDf;
}
  
export function propFunc(molecule: DG.SemanticValue, propertyCol: DG.Column, host: HTMLElement, descriptions: { [colName: string]: string }) : HTMLElement {
  const { value, addColumnIcon } = prop(
    () => propertyCol.get(molecule.cell.rowIndex),
    host,
    `Calculate ${propertyCol.name} for the whole table`,
    () => {
      const df = molecule.cell.dataFrame;
      propertyCol.name = df.columns.getUnusedName(propertyCol.name);
      propertyCol.setTag(DG.TAGS.DESCRIPTION, descriptions[propertyCol.name]);
      df.columns.add(propertyCol);
    }
  );
  return ui.divH([addColumnIcon, value], {style: {'position': 'relative'}});
}