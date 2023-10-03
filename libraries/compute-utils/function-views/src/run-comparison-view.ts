/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {CARD_VIEW_TYPE, FUNCTIONS_VIEW_TYPE, SCRIPTS_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from '../../shared-utils/consts';
import {getPropViewers} from '../../shared-utils/utils';

const RUN_NAME_COL_LABEL = 'Run name' as const;
const RUN_ID_COL_LABEL = 'RunId' as const;

/**
 * View designed to compare several FuncCalls.
 * See {@link fromComparedRuns} for limitations and creation.
 * */
export class RunComparisonView extends DG.TableView {
  /**
   * Creates a new View for provided DG.FuncCall array. Created view is derived from TableView.
   * Each row in the grid corresponds to a single DG.FuncCall.
   * @param comparedRuns Array of DG.FuncCalls to compare. Assumptions: 1) All items are homogeneous. 2) All FuncCalls are completed, i.e. all outputs are filled.
   * @param options Options for proper view placement. ParentView preserves path in breadcrumbs panel. ParentCall preserves proper placing in the toolbox
   */
  static async fromComparedRuns(
    comparedRuns: DG.FuncCall[],
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
      configFunc?: DG.Func,
    },
  ) {
    const configFunc = options.configFunc ?? comparedRuns[0].func;

    const allParamViewers = [
      ...configFunc.inputs,
      ...configFunc.outputs,
    ]
      .map((prop) => getPropViewers(prop))
      .reduce((acc, config) => {
        if (!acc[config.name])
          acc[config.name] = config.config;
        else
          acc[config.name].push(...config.config);
        return acc;
      }, {} as Record<string, Record<string, string | boolean>[]>);

    const addColumnsFromProp = (configProp: DG.Property): DG.Column[] => {
      if (configProp.propertyType === DG.TYPE.DATA_FRAME) {
        const requestedViewersConfigs = allParamViewers[configProp.name];

        const viewerColumns = requestedViewersConfigs.map((config) => {
          let columnName = configProp.caption ?? configProp.name;
          const newColumn = DG.Column.fromType(DG.TYPE.DATA_FRAME, columnName, comparedRuns.length);
          newColumn.init(
            (idx: number) => comparedRuns[idx].inputs[configProp.name] ?? comparedRuns[idx].outputs[configProp.name],
          );
          const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
          newColumn.name = unusedName;
          columnName = unusedName;
          newColumn.temp[VIEWER_PATH] = config;
          comparisonDf.columns.add(newColumn);

          return newColumn;
        });

        return viewerColumns;
      } else {
        let columnName = configProp.caption ?? configProp.name;
        //@ts-ignore
        const newColumn = DG.Column.fromType(configProp.propertyType, columnName, comparedRuns.length);
        newColumn.init(
          (idx: number) => comparedRuns[idx].inputs[configProp.name] ?? comparedRuns[idx].outputs[configProp.name],
        );
        const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
        newColumn.name = unusedName;
        columnName = unusedName;
        comparisonDf.columns.add(newColumn);
        return [newColumn];
      }
    };

    const comparisonDf = DG.DataFrame.create(comparedRuns.length);
    const uniqueRunNames = [] as string[];
    comparedRuns.forEach((run) => {
      let defaultRunName = run.options['title'] ?? `${run.func.name} - ${new Date(run.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;
      let idx = 2;
      while (uniqueRunNames.includes(defaultRunName)) {
        defaultRunName = `${run.func.name} - ${new Date(run.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})} - ${idx}`;
        idx++;
      }
      uniqueRunNames.push(defaultRunName);
    });

    comparisonDf.columns.add(DG.Column.fromStrings(
      RUN_NAME_COL_LABEL,
      uniqueRunNames,
    ));
    comparisonDf.name = options.parentCall?.func.name ? `${options.parentCall?.func.name} - comparison` : `${comparedRuns[0].func.name} - comparison`;

    configFunc.inputs.forEach((prop) => addColumnsFromProp(prop));
    configFunc.outputs.forEach((prop) => addColumnsFromProp(prop));

    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12878
    const cardView = [...grok.shell.views].find((view) => view.type === CARD_VIEW_TYPE || view.type === SCRIPTS_VIEW_TYPE || view.type === FUNCTIONS_VIEW_TYPE);
    if (cardView) grok.shell.v = cardView;

    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12879
    const tempView = grok.shell.addTableView(comparisonDf);
    tempView.temp = {'isComparison': true};
    tempView.close();

    const view = new this(tempView.dart, options);
    const comparatorFunc: string = configFunc.options['comparator'];

    setTimeout(async () => {
      view.defaultCustomize();
      if (comparatorFunc) await grok.functions.call(comparatorFunc, {'comparisonView': view});
    }, 0);

    return view;
  }

  /**
   * Constructor requires Dart object, so practically cannot be used. Use {@link fromComparedRuns} instead.
   */
  private constructor(
    dartForParent: any,
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
    } = {},
  ) {
    super(dartForParent);

    if (options.parentView) this.parentView = options.parentView;
    this.parentCall = options.parentCall || grok.functions.getCurrentCall();
  }

  private defaultCustomize(): void {
    // By default, view's content has fixed sizes - avoiding this.
    // TO DO: Find a sample for bug report
    this.root.style.removeProperty('width');
    (this.root.firstChild as HTMLElement).style.removeProperty('width');
    (this.root.firstChild!.firstChild as HTMLElement).style.removeProperty('width');

    this.grid.props.rowHeight = 180;
    this.grid.props.showAddNewRowIcon = false;
    this.grid.props.allowEdit = false;

    this.grid.columns.byName(RUN_NAME_COL_LABEL)!.width = 70;

    for (let i = 0; i < this.grid.columns.length; i++) {
      const gridCol = this.grid.columns.byIndex(i)!;
      if (gridCol.column?.temp[VIEWER_PATH]) {
        gridCol.width = 350;
        gridCol.cellType = 'html';
      };
    }

    this.grid.onCellPrepare((gc) => {
      if (gc.isColHeader || gc.isRowHeader) return;

      if (gc.tableColumn!.name === RUN_NAME_COL_LABEL)
        gc.style.textVertical = true;

      if (gc.tableColumn && gc.tableColumn.temp[VIEWER_PATH]) {
        const initialValue = gc.cell.value;

        const viewerConfig = gc.tableColumn.temp[VIEWER_PATH];
        const viewerType = viewerConfig['type'] as string;
        gc.element =
          ui.waitBox(async () => {
            const viewer = Object.values(viewerTypesMapping).includes(viewerType) ?
              DG.Viewer.fromType(viewerType, initialValue) :
              await initialValue.plot.fromType(viewerType) as DG.Viewer;

            // Workaround required since getOptions and setOptions are not symmetrical
            if (!gc.tableColumn!.temp[VIEWER_PATH]['look']) {
              viewer.setOptions(viewerConfig);
              gc.tableColumn!.temp[VIEWER_PATH] = viewer.getOptions();
            } else
              viewer.setOptions(gc.tableColumn!.temp[VIEWER_PATH]['look']);

            return viewer.root;
          });

        gc.element.style.width = '100%';
        gc.element.style.height = '100%';
      }
    });

    // Catching events to render context panel
    grok.events.onCurrentObjectChanged.subscribe(({sender}) => {
      if (
        sender instanceof DG.Column &&
        sender.type === DG.TYPE.DATA_FRAME &&
        grok.shell.tv &&
        grok.shell.tv.temp['isComparison'] &&
        [DG.VIEWER.LINE_CHART, DG.VIEWER.SCATTER_PLOT, DG.VIEWER.HISTOGRAM, DG.VIEWER.BOX_PLOT].includes(sender.temp[VIEWER_PATH]['type'])
      ) {
        grok.shell.windows.showProperties = true;

        const getAppendedDfs = (column: DG.Column) => {
          const appendedDf = column.get(0).clone() as DG.DataFrame;
          appendedDf.columns.addNew(RUN_ID_COL_LABEL, DG.TYPE.STRING).init(column.dataFrame.get(RUN_NAME_COL_LABEL, 0));

          for (let i = 1; i < column.length; i++) {
            const newRunDf = column.get(i).clone() as DG.DataFrame;
            newRunDf.columns.addNew(RUN_ID_COL_LABEL, DG.TYPE.STRING).init(column.dataFrame.get(RUN_NAME_COL_LABEL, i));

            // If one of the columns is parsed as int, it could be converted into double for proper append
            const convertibleTypes = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as DG.ColumnType[];
            for (let j = 0; j < newRunDf.columns.length; j++) {
              const newDfColumn = newRunDf.columns.byIndex(j);
              const appendDfColumn = appendedDf.columns.byIndex(j);

              if (
                newDfColumn.type !== appendDfColumn.type &&
                convertibleTypes.includes(newDfColumn.type) &&
                convertibleTypes.includes(appendDfColumn.type)
              ) {
                if (newDfColumn.type !== DG.COLUMN_TYPE.FLOAT) newRunDf.columns.replace(newDfColumn, newDfColumn.convertTo(DG.COLUMN_TYPE.FLOAT));
                if (appendDfColumn.type !== DG.COLUMN_TYPE.FLOAT) appendedDf.columns.replace(appendDfColumn, appendDfColumn.convertTo(DG.COLUMN_TYPE.FLOAT));
              }
            }

            appendedDf.append(newRunDf, true);
          }

          return appendedDf;
        };
        const unitedDf = sender.temp['unitedDf'] as DG.DataFrame ?? getAppendedDfs(sender);
        sender.temp['unitedDf'] = unitedDf;

        const config = sender.temp[VIEWER_PATH]['look'];

        // Avoiding cycling event emission
        setTimeout(() => {
          switch (sender.temp[VIEWER_PATH]['type']) {
          case DG.VIEWER.LINE_CHART:
            grok.shell.o = ui.box(unitedDf.plot.line({...config, 'split': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.SCATTER_PLOT:
            grok.shell.o = ui.box(unitedDf.plot.scatter({...config, 'color': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.HISTOGRAM:
            grok.shell.o = ui.box(unitedDf.plot.histogram({...config, 'split': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.BOX_PLOT:
            grok.shell.o = ui.box(unitedDf.plot.box({...config, 'category': RUN_ID_COL_LABEL}).root);
            break;
          }
        });
      }
    });
  }
}
