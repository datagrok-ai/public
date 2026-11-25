import * as DG from 'datagrok-api/dg';

/** Creates the viewer of the preview: of the same type and look as the [src] viewer
 *  (the inner viewer for a trellis plot), a scatter plot for a table. */
export function createPreviewViewer(src: DG.DataFrame | DG.Viewer, dataFrame: DG.DataFrame): DG.Viewer {
  const isViewer = src instanceof DG.Viewer;
  const isTrellis = isViewer && src.type === DG.VIEWER.TRELLIS_PLOT;
  const trellisInnerType = isTrellis ? src.getOptions().look['viewerType'] : null;
  const previewType = isViewer ?
    (isTrellis ? trellisInnerType : src.type) :
    DG.VIEWER.SCATTER_PLOT;
  const look = isViewer ? (isTrellis ? src.getOptions().look['innerViewerLook'] : src.getOptions().look) : null;
  const sharedOptions = {
    showDataframeFormulaLines: false,
    showViewerFormulaLines: true,
    showDataframeAnnotationRegions: false,
    showViewerAnnotationRegions: true,
    showContextMenu: false,
    axesFollowFilter: false,
    axisFont: 'normal normal 11px "Arial"',
    legendVisibility: DG.VisibilityMode.Never,
  };

  if (previewType === DG.VIEWER.LINE_CHART) {
    return DG.Viewer.lineChart(dataFrame, {
      ...sharedOptions,
      yColumnNames: look.yColumnNames?.length ? [look.yColumnNames[0]] : [],
      yAxisType: look.yAxisType,
      xAxisType: look.xAxisType,
      splitColumnNames: look.splitColumnNames,
      invertXAxis: look.invertXAxis,
      showLabels: 'Never',
      xAxisHeight: 25,
    });
  }

  if (previewType === DG.VIEWER.DENSITY_PLOT) {
    return DG.Viewer.densityPlot(dataFrame, {
      ...sharedOptions,
      xColumnName: look.xColumnName,
      yColumnName: look.yColumnName,
      xAxisType: look.xAxisType,
      yAxisType: look.yAxisType,
      invertXAxis: look.invertXAxis,
      invertYAxis: look.invertYAxis,
      showColorScale: false,
      showBinSelector: false,
    });
  }

  if (previewType === DG.VIEWER.BOX_PLOT) {
    return DG.Viewer.boxPlot(dataFrame, {
      ...sharedOptions,
      valueColumnName: look.valueColumnName,
      category1ColumnName: look.category1ColumnName,
      category2ColumnName: look.category2ColumnName,
      axisType: look.axisType,
      invertYAxis: look.invertYAxis,
      showColorSelector: false,
      showStatistics: false,
    });
  }

  if (previewType === DG.VIEWER.HISTOGRAM) {
    return DG.Viewer.histogram(dataFrame, {
      ...sharedOptions,
      valueColumnName: look.valueColumnName,
      showRangeSlider: true,
      showBinSelector: false,
      showSplitSelector: false,
      filteringEnabled: false,
      showXAxis: true,
      xAxisHeight: 25,
    });
  }

  if (previewType === DG.VIEWER.BAR_CHART) {
    return DG.Viewer.barChart(dataFrame, {
      ...sharedOptions,
      valueColumnName: look.valueColumnName,
      valueAggrType: look.valueAggrType,
      splitColumnName: look.splitColumnName,
      stackColumnName: look.stackColumnName,
      orientation: look.orientation,
      axisType: look.axisType,
      showStackSelector: false,
      // Bar chart anchors its range-slider handles at the chart edges, so with the default zero outer margins
      // the 5px handles get clipped in the small preview. Pad top/bottom so the handles fit inside the canvas.
      outerMarginTop: 15,
      outerMarginBottom: 10,
    });
  }

  return DG.Viewer.scatterPlot(dataFrame, {
    ...(look ?? {}),
    ...sharedOptions,
    showXAxis: true,
    showYAxis: true,
    showXSelector: true,
    showYSelector: true,
    yAxisType: isViewer ? look.yAxisType : 'linear',
    xAxisType: isViewer ? look.xAxisType : 'linear',
    invertXAxis: isViewer ? look.invertXAxis : false,
    invertYAxis: isViewer ? look.invertYAxis : false,
    showSizeSelector: false,
    showColorSelector: false,
    showMinMaxTickmarks: false,
    showMouseOverPoint: false,
    showCurrentPoint: false,
    zoomAndFilter: 'no action',
    xAxisHeight: 25,
  });
}
