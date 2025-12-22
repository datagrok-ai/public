import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';


async function getClickAnalysisGrid(): Promise<DG.Grid> {
  const table = await queries.getAggregatedClicks();
  table.name = 'Click Analysis';
  const descriptionCol = table.col('description');
  if (!descriptionCol)
    throw new Error('Description column is missing in the Click Analysis table');
  for (let i = 0; i < table.rowCount; i++) {
    const desc = descriptionCol.get(i);
    if (desc.includes('Shell / '))
      descriptionCol.set(i, desc.replace('Shell / ', ''));
  }

  const getPackageName = (path: string) => {
    if (!path)
      return null;
    const parts = path.split(' / ');
    let lastPackage: string | null = null;
    for (const part of parts) {
      const idx = part.indexOf(':');
      if (idx !== -1 && idx < part.length - 1 && part[idx + 1] !== ' ' && isNaN(+part[idx + 1]) &&
        !(part.includes('https:') || part.includes('http:')))
        lastPackage = part.substring(0, idx).trim();
    }
    return lastPackage;
  };

  const isPackageCol = table.columns.addNewBool('~Is Package');
  const packageNameCol = table.columns.addNewString('~Package Name');
  isPackageCol.init((i) => getPackageName(descriptionCol.get(i)) != null);
  packageNameCol.init((i) => getPackageName(descriptionCol.get(i)));

  const fillLevels = (idx: number, path: string | null) => {
    if (!path)
      return;
    const rawParts = path.split(' / ');
    const finalParts: string[] = [];

    for (const part of rawParts) {
      const colonIndex = part.indexOf(':');
      if (colonIndex !== -1 && colonIndex < part.length - 1 && part[colonIndex + 1] !== ' ')
        finalParts.push(...part.split(':'));
      else if (colonIndex !== -1 && colonIndex < part.length - 1 && part[colonIndex + 1] === ' ') {
        finalParts.push(part.substring(0, colonIndex).trim());
        finalParts.push(part.substring(colonIndex + 1).trim());
      }
      else
        finalParts.push(part);
    }

    for (let i = 0; i < finalParts.length; i++) {
      const colName = `Level ${i + 1}`;
      if (!table.columns.contains(colName))
        table.columns.addNewString(colName);
      table.col(colName)!.set(idx, finalParts[i].trim());
    }
  };

  if (table.rowCount > 0 && descriptionCol) {
    for (let i = 0; i < table.rowCount; i++)
      fillLevels(i, descriptionCol.get(i));
  }

  const grid = DG.Viewer.grid(table, {rowHeight: 20});
  grid.root.style.width = '100%';
  grid.root.style.height = '100%';

  const setGridColWidth = (colName: string, width: number) => {
    const col = grid.col(colName);
    if (col)
      col.width = width;
  };
  setGridColWidth('name', 350);
  setGridColWidth('Level 1', 120);
  setGridColWidth('Level 2', 120);
  setGridColWidth('Level 3', 120);
  grid.sort(['count'], [false]);

  return grid;
}

const filters = ui.box();
filters.style.maxWidth = '250px';
const filtersStyle = {
  columnNames: ['~Is Package', '~Package Name', 'Level 1', 'Level 2', 'Level 3', 'count'],
};


export class ClicksView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Clicks';
  }

  async initViewers(path?: string): Promise<void> {
    this.root.className = 'grok-view ui-box';
    const clickAnalysisGrid = await getClickAnalysisGrid();
    const filtersRootToAdd = DG.Viewer.filters(clickAnalysisGrid.dataFrame, filtersStyle).root;
    const filtersRootToAddChildren = Array.from(filtersRootToAdd.children);
    for (let i = 0; i < filtersRootToAddChildren.length - 1; i++)
      filtersRootToAddChildren[i].remove();
    filters.append(filtersRootToAdd);
    const treeMapViewer = DG.Viewer.treeMap(clickAnalysisGrid.dataFrame, {
      splitByColumnNames: ['Level 1', 'Level 2', 'Level 3'],
    });
    this.root.append(ui.splitH([
      filters,
      ui.splitV([
        ui.box(clickAnalysisGrid.root),
        ui.box(treeMapViewer.root),
      ]),
    ]));
  }
}
