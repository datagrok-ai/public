import {_MultiChoiceCellRenderer} from './package.g';
import {_ScatterPlotCellRenderer} from './package.g';
import {_HtmlTestCellRenderer} from './package.g';
import {_HyperlinkCellRenderer} from './package.g';
import {_BinaryImageCellRenderer} from './package.g';
import {_ImageCellRenderer} from './package.g';
import {_StarsCellRenderer} from './package.g';
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HtmlTestCellRenderer, TestCellRenderer} from './cell-types/test-cell-renderer';
import {BarCellRenderer} from './cell-types/bar-cell-renderer';

import {SparklineCellRenderer} from './sparklines/sparklines-lines';
import {BarChartCellRenderer} from './sparklines/bar-chart';
import {PieChartCellRenderer} from './sparklines/piechart';
import {RadarChartCellRender} from './sparklines/radar-chart';
import {ScatterPlotCellRenderer} from './sparklines/scatter-plot';
import {names, SparklineType, sparklineTypes, SummarySettingsBase} from './sparklines/shared';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {FormCellRenderer} from './forms/forms';
import {scWebGPUPointHitTest, scWebGPURender} from './webgpu/scatterplot';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {TagsCellRenderer} from './cell-types/tags-cell-renderer';
import {ConfidenceIntervalCellRenderer} from './cell-types/confidence-interval-cell-renderer';
export * from './package.g';
export const _package = new DG.Package();

let gpuDevice: GPUDevice | null = null;
let gpuErrorCounter = 0;

export class PackageFunctions {
  @grok.decorators.func({
    meta: {
      gridChart: 'true',
      cellType: 'bar',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static barCellRenderer() {
    return new BarCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'sparkline',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Sparklines',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static sparklineCellRenderer() {
    return new SparklineCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'barchart',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Bar Chart',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static barchartCellRenderer() {
    return new BarChartCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'piechart',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Pie Chart',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static piechartCellRenderer() {
    return new PieChartCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'radar',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Radar',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static radarCellRenderer() {
    return new RadarChartCellRender();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'smartform',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Smart Form',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static smartFormCellRenderer() {
    return new FormCellRenderer();
  }

    @grok.decorators.func({
      meta: {
        cellType: 'Tags',
        gridChart: 'true',
        virtual: 'true',
        role: 'cellRenderer'
      },
      tags: ['cellRenderer'],
      name: 'Tags',
      outputs: [{type: 'grid_cell_renderer', name: 'result'}]
    })
  static tagsCellRenderer() {
    return new TagsCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'ConfidenceInterval',
      gridChart: 'true',
      virtual: 'true',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    name: 'Confidence Interval',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
    static confidenceIntervalCellRenderer() {
      return new ConfidenceIntervalCellRenderer();
    }


  @grok.decorators.func({
    meta: {action: 'Sparklines...'},
    description: 'Adds a sparkline column for the selected columns'
  })
  static summarizeColumns(
    @grok.decorators.param({'type': 'list', 'options': {'type': 'numerical'}}) columns: DG.Column[]) {
    const table = columns[0].dataFrame;
    const name = ui.input.string('Name', {value: table.columns.getUnusedName('Summary')});
    const sparklineType = ui.input.choice('Type', {value: SparklineType.Sparkline, items: sparklineTypes});
    const columnsSelector = ui.input.columns('Columns', {value: columns, table: table,
      available: names(table.columns.numerical)});
    const hide = ui.input.bool('Hide', {value: false});
    hide.setTooltip('Hide source columns in the grid');

    function addSummaryColumn() {
      const grid = grok.shell.tv.grid;
      const left = grid.horzScroll.min;
      const columnNames = names(columnsSelector.value);
      const options = {gridColumnName: name.value, cellType: sparklineType.value!};
      const gridCol = grid.columns.add(options);
      gridCol.move(grid.columns.byName(columnNames[0])!.idx);
      gridCol.settings ??= {};
      gridCol.settings[sparklineType.value! as SparklineType] = {columnNames: columnNames};
      if (hide.value) {
        for (const name of columnNames)
          grid.columns.byName(name)!.visible = false;
      }
      grid.horzScroll.scrollTo(left);
      gridCol.scrollIntoView();
      grok.shell.o = gridCol;
    }

    DG.Dialog
      .create({title: 'Add Summary Column'})
      .add(name)
      .add(sparklineType)
      .add(columnsSelector)
      .add(hide)
      .onOK(addSummaryColumn)
      .show();
  }


  @grok.decorators.func({
    meta: {action: 'Smart form...'},
    description: 'Adds a \'form\' column for the selected columns'
  })
  static addFormColumn(
    @grok.decorators.param({'type': 'list'}) columns: DG.Column[]) {
    const table = columns[0].dataFrame;
    const name = ui.input.string('Name', {value: table.columns.getUnusedName('Form')});
    const columnsSelector = ui.input.columns('Columns', {value: columns, table: table});
    const hide = ui.input.bool('Hide', {value: false});
    hide.setTooltip('Hide source columns in the grid');

    function addSummaryColumn() {
      const grid = grok.shell.tv.grid;
      const left = grid.horzScroll.min;
      const columnNames = names(columnsSelector.value);
      const options = {gridColumnName: name.value, cellType: SparklineType.Form};
      const gridCol = grid.columns.add(options);
      gridCol.move(grid.columns.byName(columnNames[0])!.idx);
      gridCol.settings ??= {};
      gridCol.settings[SparklineType.Form] = {columnNames: columnNames};
      if (hide.value) {
        for (const name of columnNames)
          grid.columns.byName(name)!.visible = false;
      }
      grid.horzScroll.scrollTo(left);
      gridCol.scrollIntoView();
      grok.shell.o = gridCol;
    }

    DG.Dialog
      .create({title: 'Add Form'})
      .add(name)
      .add(columnsSelector)
      .add(hide)
      .onOK(addSummaryColumn)
      .show();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'testUnitsKg',
      columnTags: 'foo=bar,units=kg',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static testUnitsKgCellRenderer() {
    return new TestCellRenderer();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'testUnitsTon',
      columnTags: 'foo=bar,units=ton',
      role: 'cellRenderer'
    },
    tags: ['cellRenderer'],
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static testUnitsTonCellRenderer() {
    return new HtmlTestCellRenderer();
  }


  @grok.decorators.func()
  static demoTestUnitsCellRenderer() {
    const col1 = DG.Column.fromStrings('kg', ['a', 'b']).setTag('foo', 'bar');
    col1.meta.units = 'kg';
    col1.semType = 'test';
    const col2 = DG.Column.fromStrings('ton', ['a', 'b']).setTag('foo', 'bar');
    col2.semType = 'test';
    col2.meta.units = 'ton';
    const t = DG.DataFrame.fromColumns([col1, col2]);

    grok.shell.addTableView(t);
    grok.shell.info('Different renderers even though semantic types are the same');
  }


  @grok.decorators.autostart({tags: ['autostart']})
  static async _autoPowerGrid() {
    //PinnedUtils.registerPinnedColumns();
    DG.GridCellRenderer.register(new ScatterPlotCellRenderer());

    // handling column remove/rename in sparkline columns
    grok.events.onViewerAdded.subscribe((args) => {
      if (args.args.viewer.type !== DG.VIEWER.GRID)
        return;
      const grid = args.args.viewer as DG.Grid;
      const dataFrame = grid.dataFrame;
      if (!dataFrame)
        return;
      const getSparklineSettings =
        (gridCol: DG.GridColumn) => (gridCol.settings ?? {})[gridCol.cellType] as SummarySettingsBase;
      const findSummaryCols = (columns: (string | DG.Column)[]) => {
        const summaryCols: DG.GridColumn[] = [];
        for (let i = 1; i < grid.columns.length; i++) {
          const gridCol = grid.columns.byIndex(i)!;
          const sparklineSettings = getSparklineSettings(gridCol);
          if (sparklineTypes.includes(gridCol.cellType) && sparklineSettings?.columnNames?.length > 0 &&
            columns.some((col) => sparklineSettings.columnNames.includes(col instanceof DG.Column ? col.name : col)))
            summaryCols[summaryCols.length] = gridCol;
        }
        return summaryCols;
      };

      const colsRemovedSub = dataFrame.onColumnsRemoved.subscribe((args: DG.ColumnsArgs) => {
        const summaryCols = findSummaryCols(args.columns);
        for (const col of args.columns) {
          for (const summaryCol of summaryCols) {
            const sparklineSettings = getSparklineSettings(summaryCol);
            if (sparklineSettings.columnNames.includes(col.name))
              sparklineSettings.columnNames = sparklineSettings.columnNames.filter((name) => name !== col.name);
          }
        }
      });
      const colsRenamedSub = dataFrame.onColumnNameChanged.subscribe((args: DG.EventData) => {
        const renamedArgs: {newName: string, oldName: string} = args.args;
        const summaryCols = findSummaryCols([renamedArgs.oldName]);
        for (const summaryCol of summaryCols) {
          const sparklineSettings = getSparklineSettings(summaryCol);
          sparklineSettings.columnNames[sparklineSettings.columnNames.indexOf(renamedArgs.oldName)] =
           renamedArgs.newName;
        }
      });
      const gridDetachedSub = grid.onDetached.subscribe(() => grid.detach());
      grid.sub(colsRemovedSub);
      grid.sub(colsRenamedSub);
      grid.sub(gridDetachedSub);
    });

    if (navigator.gpu)
      gpuDevice = await getGPUDevice();
  }


  @grok.decorators.func({
    meta: {
      icon: 'files/icons/formviewer.svg',
      viewerPosition: 'bottom',
      toolbox: 'true',
      role: 'viewer'
    },
    tags: ['viewer'],
    name: 'Forms',
    description: 'Forms viewer',
    outputs: [{type: 'viewer', name: 'result'}]
  })
  static formsViewer() {
    return new FormsViewer();
  }


  @grok.decorators.panel({
    name: 'Content',
    description: 'Image content',
    meta: {role: 'widgets'},
    tags: ['widgets', 'panel']
  })
  static imgContent(
    @grok.decorators.param({'options': {'semType': 'ImageUrl'}}) imageUrl: string): DG.Widget {
    const image = new Image();
    image.src = imageUrl;
    return imageUrl !== null ? new DG.Widget(image) : new DG.Widget(ui.divText('No image available'));
  }


  @grok.decorators.func()
  static demoCellTypes() {
    const t = grok.data.demo.demog(100);
    const dis = t.col('disease')!;
    dis.set(0, 'Anxiety, Glaucoma');
    dis.set(1, 'Hepatitis A, Glaucoma');
    dis.meta.choices = ['Anxiety', 'Hepatitis A', 'Glaucoma'];
    dis.setTag(DG.TAGS.CELL_RENDERER, 'MultiChoice');

    const site = t.col('site')!;
    site.set(0, 'Buffalo, Orlando');
    site.set(1, 'Buffalo, Los Angeles');
    site.setTag(DG.TAGS.CELL_RENDERER, 'Tags');

    grok.shell.addTableView(t);
  }


  @grok.decorators.func()
  static async _scWebGPURender(
    sc: DG.ScatterPlotViewer,
    show: boolean) {
    try {
      await scWebGPURender(sc, show);
    } catch (error) {
      gpuErrorCounter++;
    }
  }


  @grok.decorators.func({
    outputs: [{type: 'int', name: 'result'}]
  })
  static async _scWebGPUPointHitTest(
    sc: DG.ScatterPlotViewer,
    pt: DG.Point) {
    let result = -1;
    try {
      result = await scWebGPUPointHitTest(sc, pt);
    } catch (error) {
      gpuErrorCounter++;
      throw error;
    }

    return result;
  }


  @grok.decorators.func()
  static isWebGPUAvailable() : boolean {
    return gpuDevice != null && !gpuErrorCounter;
  }


  @grok.decorators.func()
  static isWebGPURenderValid(sc: DG.ScatterPlotViewer) : boolean {
    return sc.props.zoomAndFilter != 'pack and zoom by filter' &&
      !sc.props.markersColumnName &&// different markers not yet supported YET
      (sc.props.jitterSize ?? 0) == 0 && (sc.props.jitterSizeY ?? 0) == 0; // jittering not supported yet
  }
}
export {_ImageCellRenderer};
export {_BinaryImageCellRenderer};
export {_HyperlinkCellRenderer};
export {_HtmlTestCellRenderer};
export {_ScatterPlotCellRenderer};
export {_MultiChoiceCellRenderer};
export {_StarsCellRenderer};
