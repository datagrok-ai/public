/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from './plate';
import {FIT_FUNCTION_4PL_REGRESSION, IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {AnalysisMappingPanel, AnalysisRequiredFields} from '../plates/views/components/analysis-mapping/analysis-mapping-panel';

export class PlateDoseRatioAnalysis {
  private static REQUIRED_FIELDS: AnalysisRequiredFields[] = [
    {
      name: 'Agonist_Concentration_M',
      required: true,
      description: 'Agonist concentration values in molar units'
    },
    {
      name: 'Antagonist_Concentration_M',
      required: true,
      description: 'Antagonist concentration values in molar units'
    },
    {
      name: 'Percent_Inhibition',
      required: true,
      description: 'Response values as percent inhibition'
    },
    {
      name: 'Agonist_ID',
      required: false,
      description: 'Optional: Identifier for the agonist compound'
    },
    {
      name: 'Antagonist_ID',
      required: false,
      description: 'Optional: Identifier for the antagonist compound'
    }
  ];

  // In PlateDoseRatioAnalysis
  static createAnalysisView(
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (target: string, source: string) => void,
    onUndo: (target: string) => void
  ): HTMLElement {
    const container = ui.divV([], 'dose-ratio-analysis-container');

    // Check if all required fields are mapped
    const requiredFields = ['Agonist_Concentration_M', 'Antagonist_Concentration_M', 'Percent_Inhibition'];
    const allRequiredMapped = requiredFields.every((field) => currentMappings.has(field));

    if (allRequiredMapped) {
      console.log('[DEBUG] All required fields mapped, creating grid...');
      const doseRatioGrid = this.createDoseRatioGrid(plate, currentMappings);

      if (doseRatioGrid) {
        console.log('[DEBUG] Grid dimensions:', doseRatioGrid.style.width, doseRatioGrid.style.height);
        console.log('[DEBUG] Grid element:', doseRatioGrid);
        console.log('[DEBUG] Container before adding grid:', container);

        container.appendChild(doseRatioGrid);

        console.log('[DEBUG] Container after adding grid:', container);
        console.log('[DEBUG] Container children count:', container.children.length);
        console.log('[DEBUG] Container dimensions:', container.style.width, container.style.height);

        // Force container to be visible
        container.style.minHeight = '400px';
        // container.style.border = '2px solid red'; // Temporary debug border
        container.style.width = '100%'; // Temporary debug border
        doseRatioGrid.style.minHeight = '400px';

        console.log('[DEBUG] Grid added to container with forced dimensions');
      }
    } else {
    // Show the mapping panel
      const mappingPanel = new AnalysisMappingPanel({
        analysisName: 'Dose Ratio',
        requiredFields: [
          {name: 'Agonist_Concentration_M', required: true},
          {name: 'Antagonist_Concentration_M', required: true},
          {name: 'Percent_Inhibition', required: true},
          {name: 'Agonist_ID', required: false},
          {name: 'Antagonist_ID', required: false}
        ],
        sourceColumns: plate.data.columns.names(),
        currentMappings,
        onMap,
        onUndo
      });

      container.appendChild(mappingPanel.getRoot());
    }

    return container;
  }


  static updateAnalysisView(
    container: HTMLElement,
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (sourceColumn: string, targetProperty: string) => void,
    onUndo: (targetProperty: string) => void
  ): void {
    const mappingPanel = (container as any).__mappingPanel as AnalysisMappingPanel;
    const resultsContainer = (container as any).__resultsContainer as HTMLElement;

    if (!mappingPanel || !resultsContainer) return;

    // Update mapping panel
    mappingPanel.updateConfig({
      sourceColumns: plate.data.columns.names(),
      currentMappings,
      onMap,
      onUndo
    });

    // Update results
    ui.empty(resultsContainer);

    const validationStatus = mappingPanel.getValidationStatus();

    if (validationStatus.isValid) {
      const chartElement = this.createDoseRatioGrid(plate, currentMappings);
      if (chartElement) {
        resultsContainer.appendChild(chartElement);
      } else {
        resultsContainer.appendChild(
          ui.divText('Unable to create dose-ratio curves with current data mapping.', 'warning-message')
        );
      }
    } else {
      const missingFieldsText = validationStatus.missingFields.join(', ');
      resultsContainer.appendChild(
        ui.divText(
          `Please map the following required fields to continue: ${missingFieldsText}`,
          'info-message'
        )
      );
    }
  }


  static createDoseRatioGrid(plate: Plate, mappings?: Map<string, string>): HTMLElement | null {
    // Use mappings to get actual column names
    console.log('[DEBUG] createDoseRatioGrid called with mappings:', mappings);

    const agonistColumn = mappings?.get('Agonist_Concentration_M') || 'Agonist_Concentration_M';
    const antagonistColumn = mappings?.get('Antagonist_Concentration_M') || 'Antagonist_Concentration_M';
    const responseColumn = mappings?.get('Percent_Inhibition') || 'Percent_Inhibition';

    console.log('[DEBUG] Looking for columns:', {agonistColumn, antagonistColumn, responseColumn});
    console.log('[DEBUG] Available columns:', plate.data.columns.names());

    // Check if required columns exist
    if (!plate.data.columns.contains(agonistColumn) ||
      !plate.data.columns.contains(antagonistColumn) ||
      !plate.data.columns.contains(responseColumn)) {
      console.warn(`[DEBUG] Missing columns. Looking for: ${agonistColumn}, ${antagonistColumn}, ${responseColumn}`);
      return null;
    }

    console.log('[DEBUG] All columns found, proceeding with grid creation...');

    const antagonistCol = plate.data.col(antagonistColumn)!;
    const agonistCol = plate.data.col(agonistColumn)!;
    const responseCol = plate.data.col(responseColumn)!;

    const uniqueConcentrations = Array.from(new Set(antagonistCol.toList()));
    const antagonistConcentrations = (uniqueConcentrations as number[]).sort((a, b) => a - b);

    const series: IFitSeries[] = [];

    for (let i = 0; i < antagonistConcentrations.length; i++) {
      const antagonistConc = antagonistConcentrations[i];
      const color = DG.Color.toHtml(DG.Color.getCategoricalColor(i));
      const currentSeries: IFitSeries = {
        name: `Antagonist: ${antagonistConc === 0 ? 'None' : `${antagonistConc.toExponential(1)} M`}`,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        points: [],
        pointColor: color,
        fitLineColor: color,
        droplines: ['IC50'],
        showPoints: 'points',
      };

      for (let j = 0; j < plate.data.rowCount; j++) {
        if (antagonistCol.get(j) === antagonistConc) {
          currentSeries.points.push({
            x: agonistCol.get(j),
            y: responseCol.get(j),
            outlier: false,
          });
        }
      }

      currentSeries.points.sort((a, b) => a.x - b.x);

      if (currentSeries.points.length > 1)
        series.push(currentSeries);
    }

    if (series.length === 0)
      return ui.divText('No data available to plot dose-ratio curves.');

    const chartData: IFitChartData = {
      chartOptions: {
        title: `Dose Ratio Analysis for ${plate.barcode}`,
        xAxisName: 'Agonist Concentration (M)',
        yAxisName: 'Percent Inhibition',
        logX: true,
        minY: -10,
        maxY: 110,
      },
      series: series,
    };

    const chartDf = DG.DataFrame.create(1);
    const curveCol = chartDf.columns.addNewString('Dose Ratio Curve');
    curveCol.init(() => JSON.stringify(chartData));
    curveCol.semType = 'fit';

    const grid = chartDf.plot.grid();

    const container = ui.div([grid.root], 'drc-grid-container');

    container.style.width = '100%';
    container.style.height = '100%';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';


    console.log('[DEBUG] Number of series created:', series.length);
    console.log('[DEBUG] Series details:', series.map((s) => ({name: s.name, pointCount: s.points.length})));

    if (series.length === 0) {
      console.log('[DEBUG] No series created - returning text message');
      return ui.divText('No data available to plot dose-ratio curves.');
    }

    // After you create the chartData, add this:
    console.log('[DEBUG] Chart data created:', chartData);

    // After you create the grid, add this:
    console.log('[DEBUG] Datagrok grid created:', grid);
    console.log('[DEBUG] Grid root:', grid.root);
    console.log('[DEBUG] Grid root style:', grid.root.style.cssText);
    ui.tools.handleResize(container, (w: number, h: number) => {
      if (w > 20 && h > 20) {
        grid.col(curveCol.name)!.width = w - 20;
        grid.props.rowHeight = h - 20;
      }
    });

    return container;
  }

  // Static method to get required fields (useful for other components)
  static getRequiredFields(): AnalysisRequiredFields[] {
    return [...this.REQUIRED_FIELDS];
  }
}
