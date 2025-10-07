/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate';
import {PlateWidget} from '../../plate-widget';
import {AbstractPlateAnalysis, IAnalysisProperty} from '../base-analysis';
import {AnalysisRequiredFields} from '../../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {FIT_FUNCTION_4PL_DOSE_RESPONSE, FIT_FUNCTION_4PL_REGRESSION, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {DrcAnalysisCoordinator} from './drc-coordinator';
import {Subscription} from 'rxjs';

export class DrcAnalysis extends AbstractPlateAnalysis {
  readonly name: string = 'DRC';
  readonly friendlyName: string = 'Dose Response';

  private coordinator?: DrcAnalysisCoordinator;
  private plateSubscription?: Subscription;

  parameters: IAnalysisProperty[] = [
    {name: 'Normalize', type: DG.TYPE.BOOL, defaultValue: true, category: 'Parameter'},
  ];
  outputs: IAnalysisProperty[] = [
    {name: 'Curve', type: DG.TYPE.STRING, category: 'Output'},
    {name: 'IC50', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'Hill Slope', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'R Squared', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'Min', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'Max', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'AUC', type: DG.TYPE.FLOAT, category: 'Output'},
  ];

  override formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame {
    console.log('--- DRC: formatResultsForGrid START ---');
    if (rawResults.rowCount === 0)
      return rawResults;

    const propsCol = rawResults.col('properties');
    if (!propsCol) {
      console.error('DRC: CRITICAL - `properties` column not found!');
      return DG.DataFrame.create(0);
    }

    const finalRows: any[] = [];
    for (let i = 0; i < rawResults.rowCount; i++) {
      const row = rawResults.row(i);
      const propsJson = propsCol.get(i);
      if (!propsJson) continue;

      try {
        const properties = JSON.parse(propsJson);
        const groupArray = row.get('group_combination');

        finalRows.push({
          'run_id': row.get('run_id'),
          'plate_id': row.get('plate_id'),
          'barcode': row.get('barcode'),
          'Compound': Array.isArray(groupArray) && groupArray.length > 0 ? groupArray[0] : null,
          'Curve': properties.Curve,
          'IC50': properties.IC50,
          'Hill Slope': properties['Hill Slope'],
          'R Squared': properties['R Squared'],
          'Min': properties.Min,
          'Max': properties.Max,
          'AUC': properties.AUC,
        });
      } catch (e) {
        console.error(`Failed to parse DRC properties JSON for row ${i}:`, propsJson, e);
      }
    }

    if (finalRows.length === 0)
      return DG.DataFrame.create(0);

    const resultDf = DG.DataFrame.fromObjects(finalRows)!;

    // Set semtypes and formats for proper rendering
    const curveCol = resultDf.col('Curve');
    if (curveCol)
      curveCol.semType = 'fit';

    const ic50Col = resultDf.col('IC50');
    if (ic50Col)
      ic50Col.meta.format = '0.00e0';

    console.log('--- DRC: formatResultsForGrid END ---');
    return resultDf;
  }


  getRequiredFields(): AnalysisRequiredFields[] {
    return [
      {name: 'Activity', required: true, description: 'Response/activity values'},
      {name: 'Concentration', required: true, description: 'Concentration/dose values'},
      {name: 'SampleID', required: true, description: 'Sample/compound identifiers'}
    ];
  }
  getSearchableProperties(): IAnalysisProperty[] {
    // DRC can search by all its outputs
    return this.outputs;
  }

  private async run(plate: Plate, params: Record<string, any>, mappings: Map<string, string>): Promise<{resultsDf: DG.DataFrame, seriesVals: Array<[string, any]>, finalValueCol: string} | null> {
    const sampleColName = mappings.get('SampleID')!;
    const concentrationColName = mappings.get('Concentration')!;
    const activityColName = mappings.get('Activity')!;

    if (!plate.data.columns.contains(sampleColName) || !plate.data.columns.contains(concentrationColName) || !plate.data.columns.contains(activityColName))
      return null;

    let finalValueCol = activityColName;
    let normed = false;
    const controlColumns = ['High Control', 'Low Control'];

    if (params['Normalize'] && controlColumns.every((c) => plate.data.col(sampleColName)!.categories.includes(c))) {
      const [lStats, hStats] = controlColumns
        .map((colName) => plate.getStatistics(activityColName, ['mean', 'std'], {match: {[sampleColName]: colName}}))
        .sort((a, b) => a.mean - b.mean);

      if (hStats.mean !== lStats.mean) {
        finalValueCol = plate.normalize(activityColName, (v) => (hStats.mean - v) / (hStats.mean - lStats.mean) * 100).name;
        normed = true;
      }
    }

    const series = plate.doseResponseSeries({
      value: finalValueCol,
      concentration: concentrationColName,
      groupBy: sampleColName
    });

    const seriesVals = Object.entries(series);
    if (seriesVals.length === 0 || !seriesVals.some(([_, s]) => s.points.length > 1))
      return null;

    const minMax = normed ? {minY: -10, maxY: 110} : {};
    const roleCol = DG.Column.string(sampleColName, seriesVals.length).init((i) => seriesVals[i][0]);
    const curveCol = DG.Column.string('Curve', seriesVals.length);
    curveCol.semType = 'fit';

    curveCol.init((i) => {
      const currentSeries = seriesVals[i][1];
      const seriesData: IFitSeries = {
        name: currentSeries.name,
        points: currentSeries.points,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        parameters: undefined,
        clickToToggle: true,
        droplines: ['IC50'],
        showPoints: 'points',
      };

      const chartData = {
        chartOptions: {
          xAxisName: concentrationColName,
          yAxisName: finalValueCol,
          logX: true,
          title: `${seriesVals[i][0]}`,
          ...minMax,
        },
        series: [seriesData],
      };

      return JSON.stringify(chartData);
    });

    const resultsDf = DG.DataFrame.fromColumns([roleCol, curveCol]);

    const statsToAdd: Record<string, string> = {'interceptX': 'IC50', 'slope': 'Hill Slope', 'rSquared': 'R Squared', 'bottom': 'Min', 'top': 'Max', 'auc': 'AUC'};

    for (const [statName, colName] of Object.entries(statsToAdd)) {
      if (this.outputs.some((o) => o.name === colName)) {
        const funcParams = {table: resultsDf, colName: 'Curve', propName: statName, seriesNumber: 0};
        const col = (DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(funcParams).callSync({processed: false})).getOutputParamValue();
        col.name = colName;
      }
    }
    if (resultsDf.col('IC50'))
        resultsDf.col('IC50')!.meta.format = 'scientific';

    return {resultsDf, seriesVals, finalValueCol};
  }

  createView(plate: Plate, plateWidget: PlateWidget, currentMappings: Map<string, string>, onMap: (t: string, s: string) => void, onUndo: (t: string) => void, onRerender?:()=>void): HTMLElement {
    return this._createStandardMappingView(plate, currentMappings, onMap, onUndo,
      (mappedPlate, mappedMappings) => {
        this.coordinator?.destroy();
        this.plateSubscription?.unsubscribe();

        const container = ui.divV([ui.loader()], 'drc-analysis-container');
        container.style.width = '100%';
        container.style.height = '100%';

        const params = {'Normalize': this.parameters.find((p) => p.name === 'Normalize')?.defaultValue ?? true};

        this.run(mappedPlate, params, mappedMappings).then((runOutput) => {
          ui.empty(container);

          if (!runOutput) {
            container.appendChild(ui.divText('Analysis failed: Not enough data to fit curves. Check mappings.', 'warning-message'));
            return;
          }

          const {resultsDf, seriesVals, finalValueCol} = runOutput;
          const resultsGrid = resultsDf.plot.grid();
          resultsGrid.props.rowHeight = 200;
          resultsGrid.root.style.width = '100%';
          resultsGrid.root.style.flexGrow = '1';

          setTimeout(() => { if (resultsGrid.col('Curve')) resultsGrid.col('Curve')!.width = 400; }, 200);

          this.coordinator = new DrcAnalysisCoordinator(
            mappedPlate, plateWidget, resultsGrid, resultsDf.col('Curve')!, seriesVals,
            {
              roleName: mappedMappings.get('SampleID')!,
              concentrationName: mappedMappings.get('Concentration')!,
              valueName: finalValueCol,
            }
          );

          this.plateSubscription = mappedPlate.onOutlierChanged.subscribe(() => {
            this.coordinator?.regenerateCurves();
          });

          const saveButton = ui.button('SAVE RESULTS', async () => await this.saveResults(mappedPlate, resultsDf, params, mappedMappings));

          saveButton.style.marginTop = '8px';

          const resultsContainer = ui.divV([
            resultsGrid.root,
            ui.div([saveButton], {style: {display: 'flex', justifyContent: 'flex-end', paddingRight: '4px'}}),
          ]);
          resultsContainer.style.cssText = 'display: flex; flex-direction: column; width: 100%; height: 100%;';
          resultsGrid.root.style.flexGrow = '1';

          container.appendChild(resultsContainer);
        });

        return container;
      }
    );
  }

  protected _getGroups(resultsDf: DG.DataFrame): { groupColumn: string, groups: string[] } {
    const groupColumn = resultsDf.columns.byIndex(0).name;
    return {groupColumn: groupColumn, groups: resultsDf.col(groupColumn)!.categories};
  }
}
