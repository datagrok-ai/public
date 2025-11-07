/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptideUtils} from '../peptideUtils';
import wu from 'wu';
import {MutationCliffs} from '../utils/types';
import {ParallelMutationCliffs} from '../utils/parallel-mutation-cliffs';
import $ from 'cash-dom';
import {PeptidesModel} from '../model';
import {extractColInfo} from '../utils/misc';
import {Subscription} from 'rxjs';
export type MutationCliffsWithMonomers = {
  cliffs: MutationCliffs,
  monomers: string[]
}

export class MutationCliffsViewer extends DG.JsViewer {
  public sequenceColumnName: string;
  public seriesColumnName: string;
  public activityColumnName: string;
  public position = 1;
  public yAxisType: 'Linear' | 'Logarithmic' = 'Linear';
  constructor() {
    super();
    this.sequenceColumnName = this.column('sequence', {semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.seriesColumnName = this.column('series', {columnTypeFilter: 'categorical', nullable: false});
    this.activityColumnName = this.column('activity', {columnTypeFilter: 'numerical', nullable: false});
    this.position = this.int('position', 1, {nullable: false, showSlider: false, min: 1, max: 100, showPlusMinus: true, description: 'Position in the sequence to analyze (1 Based).', category: 'Data'});
    this.yAxisType = this.string('yAxisType', 'Linear', {choices: ['Linear', 'Logarithmic'], description: 'Y-Axis scale type.', nullable: false, category: 'Data'}) as 'Linear' | 'Logarithmic';
  }

  onTableAttached(): void {
    super.onTableAttached();
    const seqCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    if (seqCol)
        this.getProperty('sequenceColumnName')!.set(this, seqCol.name);
    this.getProperty('activityColumnName')!.set(this, wu(this.dataFrame.columns.numerical).next().value.name);
    const seriesCol = wu(this.dataFrame.columns.categorical)?.filter((c) => c.name !== this.sequenceColumnName && c.name?.toLowerCase().includes('series')).next()?.value;
    if (seriesCol)
        this.getProperty('seriesColumnName')!.set(this, seriesCol.name);

    this.subs.push(DG.debounce(this.dataFrame.onFilterChanged, 200).subscribe(() => {
      this.clearCache();
      this.debouncedRender();
    }));
    this.debouncedRender();


    this.subs.push(grok.events.onContextMenu.subscribe((e) => {
      if (e.causedBy && e.causedBy.target && this._lineChart?.root.contains(e.causedBy.target)) {
        e.causedBy.preventDefault();
        e.causedBy.stopPropagation();
        e.causedBy.stopImmediatePropagation();
      }
    }));
  }

  private _mutationCliffsData: Promise<MutationCliffsWithMonomers> | null = null;
  private get mutationCliffsData(): Promise<MutationCliffsWithMonomers> {
    if (this._mutationCliffsData)
      return this._mutationCliffsData;

    if (!this.activityColumnName || !this.sequenceColumnName || !this.position)
      throw new Error('Activity column or Sequence column is not set, or position is invalid');

    this._mutationCliffsData = new Promise<MutationCliffsWithMonomers>(async (resolve, reject) => {
      try {
        // const seqCol = this.dataFrame.col(this.sequenceColumnName)!;
        const activityColRawData = this.dataFrame.col(this.activityColumnName)!.getRawData();

        // const seqHelper = PeptideUtils.getSeqHelper();
        // const seqHandler = seqHelper.getSeqHandler(seqCol);
        // get single position
        // const monomersAtPosition = seqHandler.getMonomersAtPosition(this.position - 1, true);
        // const monomerCol = DG.Column.fromList('string', `Monomers at pos ${this.position}`, monomersAtPosition);
        // const monomerRawData: RawColumn = {
        //   rawData: monomerCol.getRawData(),
        //   name: monomerCol.name,
        //   cat: monomerCol.categories,
        // };

        const positionColumns = this.positionColumns;
        const monomersAtPosition = positionColumns[this.position - 1].toList() as string[];
        const monomerRawData = positionColumns.map(extractColInfo);

        const mutationCliffs = await new ParallelMutationCliffs().calc(activityColRawData, monomerRawData,
          {maxMutations: 1, minActivityDelta: 0, filter: this.dataFrame.filter.anyFalse ? new Uint32Array(this.dataFrame.filter.getBuffer().buffer) : undefined, singlePosition: {position: this.position - 1}},
        );
        resolve({cliffs: mutationCliffs, monomers: monomersAtPosition});
      } catch (error) {
        reject(error);
      }
    });
    return this._mutationCliffsData;
  }

  private _dfSubs: Subscription[] = [];

  private async calculateDf() {
    const mutationCliffs = await this.mutationCliffsData;
    const uniqueIndexes = new Set<number>();
    mutationCliffs.cliffs.forEach((positionMap) => {
      positionMap.forEach((mcData) => { // should be only one position
        Array.from(mcData.entries()).forEach(([from, toIndexes]) => {
          uniqueIndexes.add(Number(from));
          toIndexes.forEach((toIdx) => uniqueIndexes.add(Number(toIdx)));
        });
      });
    });
    const indexesArray = Array.from(uniqueIndexes.values()).sort((a, b) => a - b);
    const dataFrameIndexToViewerDfIndex = new Map<number, number>();
    indexesArray.forEach((dataFrameIdx, viewerDfIdx) => {
      dataFrameIndexToViewerDfIndex.set(dataFrameIdx, viewerDfIdx);
    });
    const monomers = indexesArray.map((i) => mutationCliffs.monomers[i]);
    const activityColRawData = this.dataFrame.col(this.activityColumnName)!.getRawData();
    const activityValues = indexesArray.map((i) => activityColRawData[i] as number);
    const monomerCol = DG.Column.fromList('string', `Position ${this.position}`, monomers);
    const activityCol = DG.Column.fromList('double', this.activityColumnName, activityValues);
    // add original indexes to be able to trace back

    const cols = [monomerCol, activityCol];
    if (this.seriesColumnName) {
      const seriesColRawData = this.dataFrame.col(this.seriesColumnName)!.getRawData();
      const seriesCats = this.dataFrame.col(this.seriesColumnName)!.categories;
      const seriesValues = indexesArray.map((i) => seriesCats[seriesColRawData[i]] as string);
      const seriesCol = DG.Column.fromList('string', this.seriesColumnName, seriesValues);
      cols.push(seriesCol);
    }
    const df = DG.DataFrame.fromColumns(cols);
    // add dummy columns corresponding to other numerical columns

    Array.from(this.dataFrame.columns.numerical)
      .filter((c) => c.name !== this.activityColumnName && c.name !== this.seriesColumnName)
      .forEach((c) => {
        df.columns.addNew(c.name, c.type);
      });

    // mouse over row handling back and forth
    this._dfSubs.push(
      df.onMouseOverRowChanged.subscribe((_) => {
        if (df.mouseOverRowIdx >= 0) {
          const originalIdx = indexesArray[df.mouseOverRowIdx];
          this.dataFrame.mouseOverRowIdx = originalIdx;
        }
      }),
    );
    this._dfSubs.push(
      this.dataFrame.onMouseOverRowChanged.subscribe((_) => {
        if (this.dataFrame.mouseOverRowIdx >= 0) {
          const viewerDfIdx = dataFrameIndexToViewerDfIndex.get(this.dataFrame.mouseOverRowIdx);
          if (viewerDfIdx != undefined)
            df.mouseOverRowIdx = viewerDfIdx;
          else
            df.mouseOverRowIdx = -1;
        }
      }),
    );

    this._dfSubs.push(
      df.onCurrentRowChanged.subscribe((_) => {
        if (df.currentRowIdx >= 0) {
          const originalIdx = indexesArray[df.currentRowIdx];
          this.dataFrame.currentRowIdx = originalIdx;
        }
      }),
    );

    let firedFromViewer = false;
    let firedFromTable = false;
    // Handle selection
    this._dfSubs.push(
      DG.debounce(df.onSelectionChanged, 100).subscribe((_) => {
        const selected = df.selection;
        if (firedFromViewer) {
          firedFromViewer = false;
          return;
        }
        firedFromTable = true;
        this.dataFrame.selection.setAll(false, false);
        if (!selected.anyTrue) {
          this.dataFrame.selection.fireChanged();
          return;
        }
        for (let rowIdx = -1; (rowIdx = selected.findNext(rowIdx, true)) !== -1;) {
          const originalIdx = indexesArray[rowIdx];
          this.dataFrame.selection.set(originalIdx, true, false);
        }
        this.dataFrame.selection.fireChanged();
      }),
    );

    this._dfSubs.push(
      DG.debounce(this.dataFrame.onSelectionChanged, 100).subscribe((_) => {
        const selected = this.dataFrame.selection;
        if (firedFromTable) {
          firedFromTable = false;
          return;
        }
        firedFromViewer = true;
        df.selection.setAll(false, false);
        if (!selected.anyTrue) {
          df.selection.fireChanged();
          return;
        }
        for (let rowIdx = -1; (rowIdx = selected.findNext(rowIdx, true)) !== -1;) {
          const viewerDfIdx = dataFrameIndexToViewerDfIndex.get(rowIdx);
          if (viewerDfIdx != undefined)
            df.selection.set(viewerDfIdx, true, false);
        }
        df.selection.fireChanged();
      }),
    );

    // mouse over row GROUP handling one way
    // legend will be responsible for this
    // this._dfSubs.push(
    //   df.onMouseOverRowGroupChanged.subscribe((_) => {
    //     const func = df.rows.mouseOverRowFunc;
    //     if (!func || this.dataFrame.rowCount > 50_000)
    //       return;
    //     this.dataFrame.rows.highlight((i) => {
    //       const viewerDfIdx = dataFrameIndexToViewerDfIndex.get(i);
    //       return viewerDfIdx != undefined && func(viewerDfIdx);
    //     });
    //   }),
    // );
    return df;
  }

  private _positionColumns: DG.Column[] | null = null;
  public get positionColumns(): DG.Column[] {
    if (this._positionColumns)
      return this._positionColumns;
    // try to find model and its position columns if the SAR was run
    const peptidesModel = PeptidesModel.getInstance(this.dataFrame);
    const posCols = peptidesModel?.positionColumns;
    if (posCols && posCols.length > 0) {
      this._positionColumns = posCols;
      return this._positionColumns;
    }
    // fallback: generate columns
    const seqCol = this.dataFrame.col(this.sequenceColumnName)!;
    const seqHelper = PeptideUtils.getSeqHelper();
    const seqHandler = seqHelper.getSeqHandler(seqCol);
    const length = seqHandler.maxLength;
    const cols: DG.Column[] = [];
    for (let i = 0; i < length; i++) {
      const monomersAtPosition = seqHandler.getMonomersAtPosition(i, true);
      const monomerCol = DG.Column.fromList('string', `Position ${i + 1}`, monomersAtPosition);
      cols.push(monomerCol);
    }
    this._positionColumns = cols;
    return this._positionColumns;
  }

  private _innerDf: Promise<DG.DataFrame> | null = null;
  public get innerDf(): Promise<DG.DataFrame> {
    if (this._innerDf)
      return Promise.resolve(this._innerDf);
    this._innerDf = this.calculateDf();
    return this._innerDf;
  }

  private _lineChart: DG.LineChartViewer | null = null;

  private async render() {
    $(this.root).empty();
    if (!this.dataFrame || !this.activityColumnName || !this.sequenceColumnName || !this.position) {
      this.root.appendChild(ui.divText('Please set Activity column, Sequence column and Position properties.'));
      return;
    }

    if (this._lineChart) {
      try {
        this._lineChart.detach();
      } catch (e) {
        console.error('Error detaching previous line chart:', e);
      }
      this._lineChart = null;
    }

    ui.setUpdateIndicator(this.root, true);
    this.root.style.display = 'flex';
    this.root.style.flexDirection = 'column';

    const df = await this.innerDf;

    this._lineChart = df.plot.line({
      xColumnName: `Position ${this.position}`,
      yColumnNames: [this.activityColumnName],
      splitColumnNames: this.seriesColumnName ? [this.seriesColumnName] : [],
      legendVisibility: this.seriesColumnName ? 'Always' : 'Never',
      legendPosition: 'Right',
      showXSelector: false,
      showYSelector: true,
      showSplitSelector: false,
      xAxisLabelOrientation: 'Auto',
      axisFont: 'normal normal 14px "Roboto"',
      controlsFont: 'normal normal 14px "Roboto"',
    } as Partial<DG.ILineChartSettings>) as DG.LineChartViewer;


    this._lineChart.sub(this._lineChart.onPropertyValueChanged.subscribe((_e) => {
      if (this._lineChart?.props?.yColumnNames && this._lineChart?.props?.yColumnNames?.[0] !== this.activityColumnName) {
        const value = this._lineChart?.props?.yColumnNames?.[0];
        setTimeout(() => this.getProperty('activityColumnName')!.set(this, value), 1);
      }
    }));

    ui.setUpdateIndicator(this.root, false);
    this.root.appendChild(this._lineChart.root);

    const maxPosition = this.positionColumns.length + 1;
    const positions = new Array<number>(maxPosition - 1).fill(0).map((_, i) => i + 1).map((v) => v.toString());
    const positionInput = ui.input.choice('Position', {value: this.position.toString(), items: positions, nullable: false,
      onValueChanged: (v) => {
        if (!v)
          return;
        const intV = parseInt(v);
        this.getProperty('position')!.set(this, intV);
      },
    });

    positionInput.input.style.width = '40px';
    this.root.appendChild(ui.divH([positionInput.root], {style: {justifyContent: 'center', marginTop: '4px', width: '100%', font: 'normal normal 14px "Roboto"'}}));

    const seriesColSelector = ui.input.column('Series', {table: this.dataFrame,
      filter: (col: DG.Column) => {
        return col.isCategorical && col.name !== this.sequenceColumnName;
      }, tooltipText: 'Select column for series splitting.',
      onValueChanged: (col: DG.Column | null) => {
        const colName = col ? col.name : undefined;
        this.getProperty('seriesColumnName')!.set(this, colName);
      }, value: this.seriesColumnName ? this.dataFrame.col(this.seriesColumnName)! : undefined,
    });
    this.root.prepend(ui.divH([seriesColSelector.root], {style: {justifyContent: 'flex-end', paddingBottom: '4px', padding: '8px', width: '100%', font: 'normal normal 14px "Roboto"'}}));
  }

  private _debounceTimer: any = null;
  public debouncedRender() {
    ui.setUpdateIndicator(this.root, true);
    if (this._debounceTimer)
      clearTimeout(this._debounceTimer);
    this._debounceTimer = setTimeout(() => this.render(), 300);
  }

  private clearCache() {
    this._mutationCliffsData = null;
    this._innerDf = null;
    this._dfSubs.forEach((s) => s.unsubscribe());
    this._dfSubs = [];
  }

  detach(): void {
    super.detach();
    this._dfSubs.forEach((s) => s.unsubscribe());
    this._dfSubs = [];
  }


  onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);
    if (property?.name === 'activityColumnName' || property?.name === 'sequenceColumnName' || property?.name === 'position' || property?.name === 'seriesColumnName') {
      this.clearCache();
      this.debouncedRender();
    } if (property?.name === 'yAxisType') {
      if (this._lineChart)
        this._lineChart.props.yAxisType = this.yAxisType.toLowerCase() as DG.AxisType;
    }
  }
}
