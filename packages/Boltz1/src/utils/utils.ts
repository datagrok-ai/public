import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

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
  return ui.divH([addColumnIcon, propertyCol.get(idx)], {style: {'position': 'relative'}});
}