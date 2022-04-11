import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {PeptidesModel} from './model';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {SubstViewer} from './viewers/subst-viewer';
import {ChemPalette} from './utils/chem-palette';
import {Observable} from 'rxjs';
import {MonomerLibrary} from './monomer-library';
import {_package} from './package';
import {setAARRenderer} from './utils/cell-renderer';
import * as C from './utils/constants';
import { PeptideSpaceViewer } from './viewers/peptide-space-viewer';

type viewerTypes = SARViewer | SARViewerVertical | SubstViewer;
export class PeptidesController {
  private static _controllerName: string = 'peptidesController';
  private helpUrl = '/help/domains/bio/peptides.md';

  private _model: PeptidesModel;

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

  static get controllerName() { return PeptidesController._controllerName; }

  get dataFrame() { return this._model.dataFrame; }

  static setAARRenderer(col: DG.Column, grid: DG.Grid, grouping?: boolean) {
    return setAARRenderer(col, grid, grouping);
  }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> { return this._model.onStatsDataFrameChanged; }

  get onSARGridChanged(): Observable<DG.Grid> { return this._model.onSARGridChanged; }

  get onSARVGridChanged(): Observable<DG.Grid> { return this._model.onSARVGridChanged; }

  get onGroupMappingChanged(): Observable<StringDictionary> { return this._model.onGroupMappingChanged; }

  get onSubstTableChanged(): Observable<DG.DataFrame> { return this._model.onSubstTableChanged; }

  async updateDefault() { await this._model.updateDefault(); }

  async updateData(
    activityScaling?: string, sourceGrid?: DG.Grid, twoColorMode?: boolean, initialBitset?: DG.BitSet,
    grouping?: boolean, activityLimit?: number, maxSubstitutions?: number, isSubstitutionOn?: boolean,
  ) {
    await this._model.updateData(
      activityScaling, sourceGrid, twoColorMode, initialBitset, grouping, activityLimit, maxSubstitutions,
      isSubstitutionOn);
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

  get originalActivityColumnName(): string { return this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY]; }

  get substTooltipData() { return this._model.substTooltipData; }

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

  static get chemPalette() {
    return ChemPalette;
  }

  /**
   * Class initializer
   *
   * @param {DG.Grid} tableGrid Working talbe grid.
   * @param {DG.TableView} view Working view.
   * @param {DG.DataFrame} currentDf Working table.
   * @param {StringDictionary} options SAR viewer options
   * @param {DG.Column} col Aligned sequences column.
   * @memberof Peptides
   */
  async init(tableGrid: DG.Grid, view: DG.TableView, options: StringDictionary) {
    this.dataFrame.temp[C.EMBEDDING_STATUS] = false;
    function adjustCellSize(grid: DG.Grid) {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    }

    for (let i = 0; i < tableGrid.columns.length; i++) {
      const aarCol = tableGrid.columns.byIndex(i);
      if (aarCol && aarCol.name && aarCol.column?.semType !== C.SEM_TYPES.AMINO_ACIDS &&
        aarCol.name !== this.dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]
      )
        tableGrid.columns.byIndex(i)!.visible = false;
    }

    // const originalDfName = this.dataFrame.name;
    const dockManager = view.dockManager;

    options.title = 'Sutructure-Activity Relationship';
    const sarViewer = await this.dataFrame.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    sarViewer.helpUrl = this.helpUrl;

    const sarViewerVerticalOptions = {title: 'Most Potent Residues'};
    const sarViewerVertical =
      await this.dataFrame.plot.fromType('peptide-sar-viewer-vertical', sarViewerVerticalOptions) as SARViewerVertical;
    sarViewerVertical.helpUrl = this.helpUrl;

    const sarViewersGroup: viewerTypes[] = [sarViewer, sarViewerVertical];

    const peptideSpaceViewerOptions = {method: 't-SNE', measure: 'Levenshtein', cyclesCount: 100};
    const peptideSpaceViewer =
      await this.dataFrame.plot.fromType('peptide-space-viewer', peptideSpaceViewerOptions) as PeptideSpaceViewer;
    dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Peptide Space Viewer');

    // await createPeptideSimilaritySpaceViewer(this.dataFrame, 't-SNE', 'Levenshtein', 100, view);
    // dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Peptide Space viewer');

    dockViewers(sarViewersGroup, DG.DOCK_TYPE.RIGHT, dockManager, DG.DOCK_TYPE.DOWN);

    // const substViewerOptions = {
    //   'activityColumnName': C.COLUMNS_NAMES.ACTIVITY_SCALED,
    //   'title': 'Substitution analysis',
    // };
    // const substViewer = await this.dataFrame.plot.fromType(
    //   'substitution-analysis-viewer', substViewerOptions) as SubstViewer;
    // const substViewersGroup = [substViewer];

    tableGrid.props.allowEdit = false;
    adjustCellSize(tableGrid);

    // let isSA = false;
    // const switchViewers = ui.iconFA('toggle-on', () => {
    //   $(switchViewers).toggleClass('fa-toggle-off').toggleClass('fa-toggle-on');
    //   nodeList?.forEach((node) => {
    //     view.dockManager.close(node);
    //     node.container.destroy();
    //   });
    //   const getCurrentViewerGroup = () => isSA ? substViewersGroup : sarViewersGroup;
    //   getCurrentViewerGroup().forEach((v) => v.removeFromView());
    //   isSA = !isSA;
    //   nodeList = dockViewers(getCurrentViewerGroup(), DG.DOCK_TYPE.LEFT, dockManager, DG.DOCK_TYPE.DOWN);
    // }, 'Toggle viewer group');

    // const ribbonPanels = view.getRibbonPanels();
    // view.setRibbonPanels([[switchViewers]]);
  }
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
