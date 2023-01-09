import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {NodeType, parseNewick} from '@datagrok-libraries/bio';
import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';

/** Custom UI form for hierarchical clustering */
export async function hierarchicalClusteringUI2(df: DG.DataFrame): Promise<void> {
  async function exec(dlg: DG.Dialog) {
    const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure* Viewer');
    try {
      // TODO: Get all params from input fields
      // hierarchicalClusteringUI(df, features, distance, linkage);
      throw new Error('Not implemented');
    } catch (err) {
      const errStr = errorToConsole(err);
      grok.shell.error(errStr);
    } finally {
      pi.close();
      dlg.close();
    }
  }

  // TODO: UI to select columns for distance
  const dlg: DG.Dialog = ui.dialog(
    {title: 'Hierarchical clustering'})
    .add(ui.div([])) // TODO: UI
    .addButton('Cancel', () => { dlg.close(); })
    .addButton('Ok', async () => { await exec(dlg); })
    .show();
}

/** Creates table view with injected tree of newick result */
export async function hierarchicalClusteringUI(
  df: DG.DataFrame, colNameList: string[], distance: string, linkage: string
): Promise<void> {
  const colNameSet: Set<string> = new Set(colNameList);
  const filteredDf: DG.DataFrame = hierarchicalClusteringFilterDfForNulls(df, colNameSet);

  let tv: DG.TableView = grok.shell.getTableView(df.name);
  if (filteredDf.rowCount != df.rowCount) {
    grok.shell.warning('Hierarchical clustering analysis on data filtered out for nulls.');
    tv = grok.shell.addTableView(filteredDf);
  }
  // TODO: Filter rows with nulls in selected columns
  const preparedDf = DG.DataFrame.fromColumns(
    filteredDf.columns.toList().filter((col) => colNameSet.has(col.name)));

  const newickStr: string = await hierarchicalClusteringExec(preparedDf, distance, linkage);
  const newickRoot: NodeType = parseNewick(newickStr);
  // Fix branch_length for root node as required for hierarchical clustering result
  newickRoot.branch_length = 0;

  // empty clusterDf to stub injectTreeForGridUI2
  const clusterDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'cluster', [])]);
  injectTreeForGridUI2(tv.grid, newickRoot, undefined, 300);
}

export function hierarchicalClusteringFilterDfForNulls(df: DG.DataFrame, colNameSet: Set<string>): DG.DataFrame {
  // filteredNullsDf to open new table view
  const colList: DG.Column[] = df.columns.toList()
    .filter((col) => colNameSet.has(col.name));
  const filteredDf: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI: number) => {
    // TODO: Check nulls in columns of colNameList
    return colList.every((col) => !col.isNone(rowI));
  }));

  return filteredDf;
}

/** Runs script and returns newick result
 * @param {DG.DataFrame} preparedDf Data frame with features columns only
 */
export async function hierarchicalClusteringExec(
  preparedDf: DG.DataFrame, distance: string, linkage: string
): Promise<string> {
  // const newick: string = await hierarchicalClustering(selectedColumnsDf);
  // getNewick - is script scripts/hierarchicalClustering.py
  let newick = '';
  try {
    // Dendrogram:hierarchicalClusteringScript is the script at scripts/hierarchicalClustering.py
    newick = await grok.functions.call('Dendrogram:hierarchicalClusteringScript', {
      data: preparedDf,
      distance_name: distance,
      linkage_name: linkage,
    });
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error('Dendrogram: hierarchicalClusteringExec() ' + `${errStr}`);
    throw err;
  }

  return newick;

  // // TODO: stub
  // return new Promise<string>((resolve, reject) => {
  //   const res: string = df.getTag(thTAGS.DF_NEWICK)!;
  //   resolve(res);
  // });
}
