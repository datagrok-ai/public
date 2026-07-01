import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getGroups} from '../analysis/experiment-setup';
import {ensureDirectionColumn} from './volcano';
import {
  getIntensityColumns,
  computeMA,
  computeLoessTrend,
  computeCV,
  createMissingnessMatrix,
  unpivotIntensities,
  computeMissingBarData,
} from './qc-computations';

/** Columns added by QC computations that should be cleaned up on re-run.
 * `direction` is NOT included here: it is shared with the volcano plot, which
 * updates it in-place, and ensureDirectionColumn(...) below overwrites cells
 * rather than re-adding the column. Removing it would break any open volcano
 * viewer bound to the column reference. */
const QC_COLUMNS = ['M', 'A', 'MA_trend', 'CV_group1', 'meanInt_group1', 'CV_group2', 'meanInt_group2'];

/** Opens a QC dashboard with all QC viewers docked in a tiled layout.
 *  Requires group annotations (shows warning if missing). */
export function openQcDashboard(df: DG.DataFrame): void {
  const pi = DG.TaskBarProgressIndicator.create('Creating QC dashboard...');
  try {
    // 1. Prerequisite check: groups must be annotated
    const groups = getGroups(df);
    if (!groups) {
      grok.shell.warning('Annotate experimental groups first (Proteomics | Annotate Experiment)');
      return;
    }

    // 2. Clean up previous QC columns (allows re-running)
    for (const name of QC_COLUMNS) {
      if (df.columns.contains(name))
        df.columns.remove(name);
    }

    // 3. Get intensity columns
    const intensityCols = getIntensityColumns(df);
    if (intensityCols.length === 0) {
      grok.shell.warning('No log2 intensity columns found');
      return;
    }

    // 4. Compute MA values
    computeMA(df, groups);

    // 5. Compute loess/moving-average trend
    computeLoessTrend(df);

    // 6. Compute CV per group
    computeCV(df, groups.group1.columns, 'CV_group1', 'meanInt_group1');
    computeCV(df, groups.group2.columns, 'CV_group2', 'meanInt_group2');

    // Get the current table view
    const tv = grok.shell.tv;
    if (!tv) {
      grok.shell.warning('No table view open');
      return;
    }

    // 7. MA scatter plot
    const maPlotOptions: Record<string, any> = {
      xColumnName: 'A',
      yColumnName: 'M',
      title: 'MA Plot',
      markerDefaultSize: 2,
    };

    // Conditional coloring by direction when DE is available
    if (df.getTag('proteomics.de_complete') === 'true') {
      try {
        const dirColName = ensureDirectionColumn(df, 1.0, 0.05);
        maPlotOptions['colorColumnName'] = dirColName;
      } catch (_e) {
        // DE columns may not exist despite tag -- skip coloring
      }
    }

    const maPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, maPlotOptions);

    // Add M=0 reference line for MA plot. Drop any prior M=0 lines first so
    // repeated QC dashboard opens don't stack duplicates on the DataFrame.
    const maFormula = '${M} = 0';
    df.meta.formulaLines.items = df.meta.formulaLines.items.filter((line) => line.formula !== maFormula);
    df.meta.formulaLines.addLine({
      formula: maFormula,
      color: '#888888',
      width: 1,
      visible: true,
    });

    // 7b. MA trend line as a separate scatter plot with line markers
    const maTrend = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {
      xColumnName: 'A',
      yColumnName: 'MA_trend',
      title: 'MA Trend',
      markerDefaultSize: 1,
    });

    // 8. CV scatter plot (group1 by default)
    const cvPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {
      xColumnName: 'meanInt_group1',
      yColumnName: 'CV_group1',
      title: 'CV Plot',
      markerDefaultSize: 2,
    });

    // 9. Sample correlation plot
    const corrPlot = tv.addViewer(DG.VIEWER.CORR_PLOT, {
      xColumnNames: intensityCols,
      yColumnNames: intensityCols,
      title: 'Sample Correlation',
    });

    // 10. Missing values heatmap (separate DataFrame -- no addTable to avoid circular JSON)
    const missDf = createMissingnessMatrix(df, intensityCols);
    const missGrid = DG.Viewer.grid(missDf);
    missGrid.props.allowColSelection = false;

    // 11. Missing values bar chart (separate DataFrame -- use factory to bind correctly)
    const barDf = computeMissingBarData(df, intensityCols, groups);
    const missBar = DG.Viewer.barChart(barDf, {
      splitColumnName: 'Sample',
      valueColumnName: 'MissingPct',
      stackColumnName: 'Group',
    } as any);

    // 12. Box plot for intensity distributions (separate DataFrame -- use factory to bind correctly)
    const longDf = unpivotIntensities(df, intensityCols);
    const boxPlot = DG.Viewer.boxPlot(longDf, {
      valueColumnName: 'Intensity',
      categoryColumnName: 'Sample',
    } as any);

    // 13. Dock layout -- arrange viewers in a tiled grid
    // Top row: MA plot + MA trend + CV plot
    // Bottom row: Correlation + Missing heatmap + Missing bar + Box plot
    const dm = tv.dockManager;
    const maNode = dm.dock(maPlot, DG.DOCK_TYPE.RIGHT, null, 'MA Plot', 0.5);
    dm.dock(maTrend, DG.DOCK_TYPE.DOWN, maNode, 'MA Trend', 0.3);
    dm.dock(cvPlot, DG.DOCK_TYPE.RIGHT, maNode, 'CV Plot', 0.5);
    const corrNode = dm.dock(corrPlot, DG.DOCK_TYPE.DOWN, null, 'Sample Correlation', 0.4);
    const missGridNode = dm.dock(missGrid, DG.DOCK_TYPE.RIGHT, corrNode, 'Missing Values', 0.5);
    const missBarNode = dm.dock(missBar, DG.DOCK_TYPE.RIGHT, missGridNode, 'Missing % per Sample', 0.5);
    dm.dock(boxPlot, DG.DOCK_TYPE.RIGHT, missBarNode, 'Intensity Distributions', 0.5);
  } finally {
    pi.close();
  }
}
