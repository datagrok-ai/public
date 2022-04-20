import * as DG from 'datagrok-api/dg';
import {PeptidesModel} from './model';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {ChemPalette} from './utils/chem-palette';
import {Observable} from 'rxjs';
import {MonomerLibrary} from './monomer-library';
import {_package} from './package';
import {setAARRenderer} from './utils/cell-renderer';
import * as C from './utils/constants';
import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';

type viewerTypes = SARViewer | SARViewerVertical;
export class PeptidesController {
  private static _controllerName: string = 'peptidesController';
  private helpUrl = '/help/domains/bio/peptides.md';

  private _model: PeptidesModel;
  sarViewer!: SARViewer;
  sarViewerVertical!: SARViewerVertical;

  private constructor(dataFrame: DG.DataFrame) {
    this._model = PeptidesModel.getInstance(dataFrame);
    // this.getOrInitModel(this._dataFrame);
  }

  static async getInstance(dataFrame: DG.DataFrame): Promise<PeptidesController> {
    dataFrame.temp[PeptidesController.controllerName] ??= new PeptidesController(dataFrame);
    if (dataFrame.temp[MonomerLibrary.id] === null) {
      const sdf = await _package.files.readAsText('HELMMonomers_June10.sdf');
      dataFrame.temp[MonomerLibrary.id] ??= new MonomerLibrary(sdf);
    }
    return dataFrame.temp[PeptidesController.controllerName];
  }

  static get controllerName() {return PeptidesController._controllerName;}

  get dataFrame() {return this._model.dataFrame;}

  static setAARRenderer(col: DG.Column, grid: DG.Grid, grouping?: boolean) {
    return setAARRenderer(col, grid, grouping);
  }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {return this._model.onStatsDataFrameChanged;}

  get onSARGridChanged(): Observable<DG.Grid> {return this._model.onSARGridChanged;}

  get onSARVGridChanged(): Observable<DG.Grid> {return this._model.onSARVGridChanged;}

  get onGroupMappingChanged(): Observable<StringDictionary> {return this._model.onGroupMappingChanged;}

  get onSubstTableChanged(): Observable<DG.DataFrame> {return this._model.onSubstTableChanged;}

  async updateDefault() {await this._model.updateDefault();}

  get sarGrid() {return this._model.sarGrid;}

  get sarVGrid() {return this._model.sarVGrid;}

  get sourceGrid() {return this._model._sourceGrid!; }

  async updateData(
    activityScaling?: string, sourceGrid?: DG.Grid, twoColorMode?: boolean, grouping?: boolean, activityLimit?: number,
    maxSubstitutions?: number, isSubstitutionOn?: boolean, filterMode?: boolean,
  ) {
    await this._model.updateData(
      activityScaling, sourceGrid, twoColorMode, grouping, activityLimit, maxSubstitutions, isSubstitutionOn,
      filterMode);
  }

  static async scaleActivity(
    activityScaling: string, df: DG.DataFrame, originalActivityName?: string,
  ): Promise<[DG.DataFrame, string]> {
    // const df = sourceGrid.dataFrame!;
    let currentActivityColName = originalActivityName ?? C.COLUMNS_NAMES.ACTIVITY;
    const flag = (df.columns as DG.ColumnList).names().includes(currentActivityColName) &&
      currentActivityColName === originalActivityName;
    currentActivityColName = flag ? currentActivityColName : C.COLUMNS_NAMES.ACTIVITY;
    const tempDf = df.clone(null, [currentActivityColName]);

    let formula = '${' + currentActivityColName + '}';
    let newColName = originalActivityName ?? df.temp[C.COLUMNS_NAMES.ACTIVITY] ?? currentActivityColName;
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

    await (tempDf.columns as DG.ColumnList).addNewCalculated(C.COLUMNS_NAMES.ACTIVITY_SCALED, formula);

    return [tempDf, newColName];
  }

  get originalActivityColumnName(): string {return this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY];}

  get substTooltipData() {return this._model.substTooltipData;}

  static splitAlignedPeptides(peptideColumn: DG.Column, filter: boolean = true): [DG.DataFrame, number[]] {
    const splitPeptidesArray: string[][] = [];
    let currentSplitPeptide: string[];
    let modeMonomerCount = 0;
    let currentLength;
    const colLength = peptideColumn.length;

    // splitting data
    const monomerLengths: {[index: string]: number} = {};
    for (let i = 0; i < colLength; i++) {
      currentSplitPeptide = peptideColumn.get(i).split('-').map((value: string) => value ? value : '-');
      splitPeptidesArray.push(currentSplitPeptide);
      currentLength = currentSplitPeptide.length;
      monomerLengths[currentLength + ''] =
        monomerLengths[currentLength + ''] ? monomerLengths[currentLength + ''] + 1 : 1;
    }
    //@ts-ignore: what I do here is converting string to number the most effective way I could find. parseInt is slow
    modeMonomerCount = 1 * Object.keys(monomerLengths).reduce((a, b) => monomerLengths[a] > monomerLengths[b] ? a : b);

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
    columnNames.splice(0, 0, 'N-terminal');
    columnNames.push('C-terminal');

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

  static get chemPalette() { return ChemPalette; }

  syncProperties(isSourceSAR = true) {
    const sourceViewer = isSourceSAR ? this.sarViewer : this.sarViewerVertical;
    const targetViewer = isSourceSAR ? this.sarViewerVertical : this.sarViewer;
    const properties = sourceViewer.props.getProperties();
    for (const property of properties)
      targetViewer.props.set(property.name, property.get(sourceViewer));
  }

  modifyOrCreateSplitCol(aar: string, position: string, notify: boolean = true) {
    this._model.modifyOrCreateSplitCol(aar, position);
    if (notify)
      this._model.fireBitsetChanged();
  }

  setSARGridCellAt(aar: string, position: string) {
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
  async init(sourceGrid: DG.Grid, currentView: DG.TableView, options: StringDictionary) {
    this.dataFrame.temp[C.EMBEDDING_STATUS] = false;
    function adjustCellSize(grid: DG.Grid) {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    }

    for (let i = 0; i < sourceGrid.columns.length; i++) {
      const aarCol = sourceGrid.columns.byIndex(i);
      if (aarCol && aarCol.name && aarCol.column?.semType !== C.SEM_TYPES.AMINO_ACIDS &&
        aarCol.name !== this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]
      )
        sourceGrid.columns.byIndex(i)!.visible = false;
    }

    await this.updateData(options.scaling, sourceGrid, false, false, 1, 2, false, false);

    const dockManager = currentView.dockManager;

    this.sarViewer = await this.dataFrame.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    this.sarViewer.helpUrl = this.helpUrl;

    this.sarViewerVertical =
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

    this._model.sarGrid.invalidate();
    this._model.sarVGrid.invalidate();
  }

  invalidateSourceGrid() { this.sourceGrid.invalidate(); }
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
