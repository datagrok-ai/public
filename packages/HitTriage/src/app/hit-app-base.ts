import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {ComputeFunctions} from './types';
import {HitTriageComputeFunctionTag} from './consts';

export class HitAppBase<T> {
  public dataFrame?: DG.DataFrame;
  public template?: T;
  public baseUrl!: string;
  public computeFunctions: Promise<ComputeFunctions>;
  protected isJoining = false;
  // public layouts: Promise<DG.ViewLayout[]>;
  constructor(public parentCall: DG.FuncCall) {
    this.resetBaseUrl();
    this.computeFunctions = new Promise<ComputeFunctions>(async (resolve) => {
      const functions = DG.Func.find({tags: [HitTriageComputeFunctionTag]})
      //TODO: remove this when everyone is upgraded to support DG.FUNC.FIND returning all functions, queries and scripts
        .filter((f) => f.type === 'function-package');
      const scripts = await grok.dapi.scripts.include('params').filter(`#${HitTriageComputeFunctionTag}`).list();
      const queries = await grok.dapi.queries.include('params,connection')
        .filter(`#${HitTriageComputeFunctionTag}`).list();
      resolve({functions, scripts, queries});
    });

    // this.layouts = new Promise<DG.ViewLayout[]>(async (resolve) => {
    //   // we need all layouts, as applicable ones might not be enaugh
    //   const layouts = await grok.dapi.layouts.list();
    //   resolve(layouts);
    // });
  }

  public resetBaseUrl() {
    const href = window.location.href;
    const urlObj = new URL(href);
    this.baseUrl = urlObj.origin + urlObj.pathname;
  }

  protected getFilterType(colName: string): DG.FILTER_TYPE {
    const col = this.dataFrame!.col(colName);
    if (col?.semType === DG.SEMTYPE.MOLECULE)
      return DG.FILTER_TYPE.SUBSTRUCTURE;
    if (col?.type === DG.COLUMN_TYPE.BOOL)
      return DG.FILTER_TYPE.BOOL_COLUMNS;
    if (col?.type === DG.COLUMN_TYPE.STRING)
      return DG.FILTER_TYPE.CATEGORICAL;
    return DG.FILTER_TYPE.HISTOGRAM;
  }

  public download(df: DG.DataFrame, name: string, onlyFiltered = false): void {
    const element = document.createElement('a');
    const result = DG.DataFrame.fromColumns(df.columns.toList().filter((c) => !c.name.startsWith('~')))
      .toCsv({filteredRowsOnly: onlyFiltered});
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', name + '.csv');
    element.click();
  }

  public async getNewCampaignName(folderName: string, templateKey: string) {
    const templateCampaigns = (await _package.files.list(folderName))
      .map((file) => file.name)
      .filter((name) => name.startsWith(templateKey));
    if (templateCampaigns.length === 0)
      return templateKey + '-1';
    const postFixes = templateCampaigns.map((c) => c.split('-')[1]).filter(Boolean).map((c) => parseInt(c, 10)).sort();
    return templateKey + '-' + ((postFixes[postFixes.length - 1] + 1).toString());
  }

  protected setBaseUrl() {
    const title = document.title;
    if (history.replaceState) {
      const obj = {
        Title: title,
        Url: this.baseUrl,
      };
      history.replaceState(obj, obj.Title, obj.Url);
    }
  }

  /// creates a union of two dataframes based on a molecule column and adds to second dataframe
  protected unionDataframes(df1: DG.DataFrame, df2: DG.DataFrame, molColName: string) {
    const df1MolCol: string[] | undefined = df1.col(molColName)?.toList()?.filter((it) => !!it)
      .map((it) => DG.chem.convert(it, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles));
    const df2MolCol: string[] | undefined = df2.col(molColName)?.toList()?.filter((it) => !!it)
      .map((it) => DG.chem.convert(it, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles));
    if (!df1MolCol || !df2MolCol)
      throw new Error('Molecule column not found');
    this.isJoining = true;
    try {
    // first check that all columns are there
      for (const col of df1.columns) {
        if (!df2.columns.contains(col.name))
          df2.columns.addNew(col.name, col.type);
      }
      for (let i = 0; i < df1MolCol.length; i++) {
        if (!df2MolCol.includes(df1MolCol[i])) {
          df2.rows.addNew(null, true);
          for (const col of df1.columns) {
            const value = col.get(i);
            df2.col(col.name)?.set(df2.rowCount - 1, value, false);
          }
        }
      }
    } finally {
      setTimeout(() => {
        this.isJoining = false;
      }, 500);
    }
    return df2;
  }
}
