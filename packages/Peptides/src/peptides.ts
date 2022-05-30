import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {PeptidesModel} from './model';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {ChemPalette} from './utils/chem-palette';
import {Observable} from 'rxjs';
import {MonomerLibrary} from './monomer-library';
import {_package} from './package';
import {setAARRenderer} from './utils/cell-renderer';
import * as C from './utils/constants';
import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';
import {FilteringStatistics} from './utils/filtering-statistics';
import * as type from './utils/types';

type viewerTypes = SARViewer | SARViewerVertical;
export class PeptidesController {
  private static _controllerName: string = 'peptidesController';
  private helpUrl = '/help/domains/bio/peptides.md';

  private _model: PeptidesModel;
  sarViewer!: SARViewer;
  sarViewerVertical!: SARViewerVertical;
  isInitialized = false;

  private constructor(dataFrame: DG.DataFrame) {
    this._model = PeptidesModel.getInstance(dataFrame);
  }

  static async getInstance(dataFrame: DG.DataFrame): Promise<PeptidesController> {
    dataFrame.temp[PeptidesController.controllerName] ??= new PeptidesController(dataFrame);
    if (dataFrame.temp[MonomerLibrary.id] === null) {
      const sdf = await _package.files.readAsText('HELMMonomers_June10.sdf');
      dataFrame.temp[MonomerLibrary.id] ??= new MonomerLibrary(sdf);
    }
    return dataFrame.temp[PeptidesController.controllerName] as PeptidesController;
  }

  static get controllerName(): string {return PeptidesController._controllerName;}

  static get chemPalette(): typeof ChemPalette {return ChemPalette;}

  get dataFrame(): DG.DataFrame {return this._model.dataFrame;}

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {return this._model.onStatsDataFrameChanged;}

  get onSARGridChanged(): Observable<DG.Grid> {return this._model.onSARGridChanged;}

  get onSARVGridChanged(): Observable<DG.Grid> {return this._model.onSARVGridChanged;}

  get onSubstTableChanged(): Observable<type.SubstitutionsInfo> {return this._model.onSubstTableChanged;}

  async updateDefault(): Promise<void> {await this._model.updateDefault();}

  get sarGrid(): DG.Grid {return this._model._sarGrid;}

  get sarVGrid(): DG.Grid {return this._model._sarVGrid;}

  get sourceGrid(): DG.Grid {return this._model._sourceGrid;}

  get originalActivityColumnName(): string {return this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY] as string;}

  // get substTooltipData() {return this._model.substTooltipData;}

  static setAARRenderer(col: DG.Column, grid: DG.Grid): void {
    return setAARRenderer(col, grid);
  }

  static async scaleActivity(
    activityScaling: string, df: DG.DataFrame, originalActivityName?: string, cloneBitset = false,
  ): Promise<[DG.DataFrame, string]> {
    let currentActivityColName = originalActivityName ?? C.COLUMNS_NAMES.ACTIVITY;
    const flag = df.columns.names().includes(currentActivityColName) &&
      currentActivityColName === originalActivityName;
    currentActivityColName = flag ? currentActivityColName : C.COLUMNS_NAMES.ACTIVITY;
    const tempDf = df.clone(cloneBitset ? df.filter : null, [currentActivityColName]);

    let formula = '${' + currentActivityColName + '}';
    let newColName = 'activity'; //originalActivityName ?? df.temp[C.COLUMNS_NAMES.ACTIVITY] ?? currentActivityColName;
    switch (activityScaling) {
    case 'none':
      break;
    case 'lg':
      formula = `Log10(${formula})`;
      newColName = `Log10(${newColName})`;
      break;
    case '-lg':
      formula = `-1*Log10(${formula})`;
      newColName = `-Log10(${newColName})`;
      break;
    default:
      throw new Error(`ScalingError: method \`${activityScaling}\` is not available.`);
    }

    await tempDf.columns.addNewCalculated(C.COLUMNS_NAMES.ACTIVITY_SCALED, formula);
    df.tags['scaling'] = activityScaling;

    return [tempDf, newColName];
  }

  static splitAlignedPeptides(peptideColumn: DG.Column, filter: boolean = true): [DG.DataFrame, number[]] {
    const separator = peptideColumn.tags[C.TAGS.SEPARATOR];
    const splitPeptidesArray: string[][] = [];
    let currentSplitPeptide: string[];
    let modeMonomerCount = 0;
    let currentLength;
    const colLength = peptideColumn.length;

    // splitting data
    const monomerLengths: {[index: string]: number} = {};
    for (let i = 0; i < colLength; i++) {
      currentSplitPeptide = peptideColumn.get(i).split(separator).map((value: string) => value ? value : '-');
      splitPeptidesArray.push(currentSplitPeptide);
      currentLength = currentSplitPeptide.length;
      monomerLengths[currentLength + ''] =
        monomerLengths[currentLength + ''] ? monomerLengths[currentLength + ''] + 1 : 1;
    }
    modeMonomerCount =
      parseInt(Object.keys(monomerLengths).reduce((a, b) => monomerLengths[a] > monomerLengths[b] ? a : b));

    // making sure all of the sequences are of the same size
    // and marking invalid sequences
    let nTerminal: string;
    const invalidIndexes: number[] = [];
    let splitColumns: string[][] = Array.from({length: modeMonomerCount}, (_) => []);
    modeMonomerCount--; // minus N-terminal
    for (let i = 0; i < colLength; i++) {
      currentSplitPeptide = splitPeptidesArray[i];
      nTerminal = currentSplitPeptide.pop()!; // it is guaranteed that there will be at least one element
      currentLength = currentSplitPeptide.length;
      if (currentLength !== modeMonomerCount)
        invalidIndexes.push(i);

      for (let j = 0; j < modeMonomerCount; j++)
        splitColumns[j].push(j < currentLength ? currentSplitPeptide[j] : '-');

      splitColumns[modeMonomerCount].push(nTerminal);
    }
    modeMonomerCount--; // minus C-terminal

    //create column names list
    const columnNames = Array.from({length: modeMonomerCount}, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1 }`);
    columnNames.splice(0, 0, 'N');
    columnNames.push('C');

    // filter out the columns with the same values
    if (filter) {
      splitColumns = splitColumns.filter((positionArray, index) => {
        const isRetained = new Set(positionArray).size > 1;
        if (!isRetained)
          columnNames.splice(index, 1);

        return isRetained;
      });
    }

    return [
      DG.DataFrame.fromColumns(splitColumns.map((positionArray, index) => {
        return DG.Column.fromList('string', columnNames[index], positionArray);
      })),
      invalidIndexes,
    ];
  }

  async updateData(
    activityScaling?: string, sourceGrid?: DG.Grid, twoColorMode?: boolean, activityLimit?: number,
    maxSubstitutions?: number, isSubstitutionOn?: boolean, filterMode?: boolean,
  ): Promise<void> {
    filterMode ??= false;
    await this._model.updateData(
      activityScaling, sourceGrid, twoColorMode, activityLimit, maxSubstitutions, isSubstitutionOn, filterMode);
  }

  getSubstitutions(): type.SubstitutionsInfo {
    return this._model.substitutionsInfo;
  }

  getCurrentAARandPos(): {aar: string, pos: string} {
    return {aar: this.dataFrame.getTag(C.TAGS.AAR) ?? 'All', pos: this.dataFrame.getTag(C.TAGS.POSITION) ?? 'All'};
  }

  assertVar(variable: string, init = false): boolean {
    //@ts-ignore
    let foundVariable: any = this[variable];
    if (!foundVariable && init) {
      //@ts-ignore
      this[variable] = foundVariable = this.dataFrame.temp[variable];
    }

    const assertionResult = foundVariable ? true : false;
    if (init && !assertionResult)
      throw new Error(`Variable assertion error: variable '${variable}' is not found in dataFrame`);

    return assertionResult;
  }

  assertVariables(variables: string[], init = false): boolean {
    let result = true;
    for (const variable of variables)
      result &&= this.assertVar(variable, init);

    return result;
  }

  syncProperties(isSourceSAR = true): void {
    this.assertVariables(['sarViewer', 'sarViewerVertical'], true);
    const sourceViewer = isSourceSAR ? this.sarViewer : this.sarViewerVertical;
    const targetViewer = isSourceSAR ? this.sarViewerVertical : this.sarViewer;
    const properties = sourceViewer.props.getProperties();
    for (const property of properties)
      targetViewer.props.set(property.name, property.get(sourceViewer));
  }

  modifyOrCreateSplitCol(aar: string, position: string, notify: boolean = true): void {
    this._model.modifyOrCreateSplitCol(aar, position);
    if (notify)
      this._model.fireBitsetChanged();
  }

  setSARGridCellAt(aar: string, position: string): void {
    const sarDf = this.sarGrid.dataFrame;
    const aarCol = sarDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    const aarColLen = aarCol.length;
    let index = -1;
    for (let i = 0; i < aarColLen; i++) {
      if (aarCol.get(i) === aar) {
        index = i;
        break;
      }
    }
    position = position === C.CATEGORIES.ALL ? C.COLUMNS_NAMES.AMINO_ACID_RESIDUE : position;
    sarDf.currentCell = sarDf.cell(index, position);
  }

  /**
   * Class initializer
   *
   * @param {DG.Grid} sourceGrid Working talbe grid.
   * @param {DG.TableView} currentView Working view.
   * @param {DG.DataFrame} currentDf Working table.
   * @param {StringDictionary} options SAR viewer options
   * @param {DG.Column} col Aligned sequences column.
   * @memberof Peptides
   */
  async init(table: DG.DataFrame): Promise<void> {
    if (this.isInitialized)
      return;
    this.isInitialized = true;
    //calculate initial stats
    const stats = new FilteringStatistics();
    const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    stats.setData(activityScaledCol.getRawData() as Float32Array);
    stats.setMask(table.selection);
    table.temp[C.STATS] = stats;

    //set up views
    let currentView = grok.shell.v as DG.TableView;
    if (currentView.dataFrame.tags['isPeptidesAnalysis'] !== 'true')
      currentView = grok.shell.addTableView(table);
    const sourceGrid = currentView.grid;
    sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED)!.name = table.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED];
    sourceGrid.columns.setOrder([table.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]]);

    this.dataFrame.temp[C.EMBEDDING_STATUS] = false;
    const adjustCellSize = (grid: DG.Grid): void => {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    };

    for (let i = 0; i < sourceGrid.columns.length; i++) {
      const aarCol = sourceGrid.columns.byIndex(i);
      if (aarCol && aarCol.name && aarCol.column?.semType !== C.SEM_TYPES.AMINO_ACIDS &&
        aarCol.name !== this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]
      )
        sourceGrid.columns.byIndex(i)!.visible = false;
    }

    const options = {scaling: table.tags['scaling']};
    await this.updateData(table.tags['scaling'], sourceGrid, false, 1, 2, false, false);

    const dockManager = currentView.dockManager;

    this.dataFrame.temp['sarViewer'] = this.sarViewer =
      await this.dataFrame.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    this.sarViewer.helpUrl = this.helpUrl;

    this.dataFrame.temp['sarViewerVertical'] = this.sarViewerVertical =
      await this.dataFrame.plot.fromType('peptide-sar-viewer-vertical', options) as SARViewerVertical;
    this.sarViewerVertical.helpUrl = this.helpUrl;

    const sarViewersGroup: viewerTypes[] = [this.sarViewer, this.sarViewerVertical];

    const peptideSpaceViewerOptions = {method: 't-SNE', measure: 'Levenshtein', cyclesCount: 100};
    const peptideSpaceViewer =
      await this.dataFrame.plot.fromType('peptide-space-viewer', peptideSpaceViewerOptions) as PeptideSpaceViewer;
    dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Peptide Space Viewer');

    dockViewers(sarViewersGroup, DG.DOCK_TYPE.RIGHT, dockManager, DG.DOCK_TYPE.DOWN);

    sourceGrid.props.allowEdit = false;
    adjustCellSize(sourceGrid);

    this._model.invalidateGrids();
  }

  invalidateSourceGrid(): void {this.sourceGrid.invalidate();}
}

function dockViewers(
  viewerList: viewerTypes[], attachDirection: DG.DockType, dockManager: DG.DockManager,
  initialAttachDirection?: DG.DockType): DG.DockNode[] | null {
  const viewerListLength = viewerList.length;
  if (viewerListLength === 0)
    return null;

  let currentViewer = viewerList[0];
  const nodeList = [dockManager.dock(currentViewer, initialAttachDirection, null, currentViewer.name ?? '')];
  const ratio = 1 / viewerListLength;

  for (let i = 1; i < viewerListLength; i++) {
    currentViewer = viewerList[i];
    nodeList.push(dockManager.dock(currentViewer, attachDirection, nodeList[i - 1], currentViewer.name ?? '', ratio));
  }
  return nodeList;
}
