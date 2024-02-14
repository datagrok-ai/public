/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {filter} from 'rxjs/operators';
import {CARD_VIEW_TYPE, FUNCTIONS_VIEW_TYPE,
  SCRIPTS_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from '../../shared-utils/consts';
import {getDfFromRuns} from './shared/utils';
import {RUN_NAME_COL_LABEL} from '../../shared-utils/consts';

const api: any = window;

/**
 * View designed to compare several FuncCalls.
 * See {@link fromComparedRuns} for limitations and creation.
 * */
export class RunComparisonView extends DG.TableView {
  /**
   * Creates a new View for provided DG.FuncCall array. Created view is derived from TableView.
   * Each row in the grid corresponds to a single DG.FuncCall.
   * @param comparedRuns Array of DG.FuncCalls to compare.
   * Assumptions: 1) All items are homogeneous. 2) All FuncCalls are completed, i.e. all outputs are filled.
   * @param options Options for proper view placement.
   * ParentView preserves path in breadcrumbs panel. ParentCall preserves proper placing in the toolbox
   */
  static async fromComparedRuns(
    comparedRuns: DG.FuncCall[],
    func: DG.Func,
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
    } = {parentView: undefined, parentCall: undefined},
  ) {
    const comparisonDf = getDfFromRuns(
      comparedRuns,
      func,
      options,
    );

    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12878
    const cardView = [...grok.shell.views].find((view) =>
      view.type === CARD_VIEW_TYPE ||
      view.type === SCRIPTS_VIEW_TYPE ||
      view.type === FUNCTIONS_VIEW_TYPE);
    if (cardView) grok.shell.v = cardView;

    const defaultView = new this(comparisonDf, options);

    return defaultView;
  }

  /**
   * Constructor requires Dart object, so practically cannot be used. Use {@link fromComparedRuns} instead.
   */
  private constructor(
    dataFrame: DG.DataFrame,
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
    } = {},
  ) {
    super(api.grok_TableView(dataFrame.dart, false));

    if (options.parentView) this.parentView = options.parentView;
    this.parentCall = options.parentCall || grok.functions.getCurrentCall();
  }

  private showAllColumns() {
    const allColumns = this.dataFrame.columns.names();
    this.grid.columns.setVisible(allColumns);
  }

  private showChangingColumns() {
    const changingColumns =
        this.dataFrame.columns.names().filter((name) =>
          (DG.TYPES_SCALAR.has(this.dataFrame.getCol(name).type as DG.TYPE) &&
          this.dataFrame.getCol(name).categories.length > 1) ||
          !DG.TYPES_SCALAR.has(this.dataFrame.getCol(name).type as DG.TYPE),
        );
    this.grid.columns.setVisible(changingColumns);
  }

  public defaultCustomize(): void {
    const showChangingColumnsIcon = ui.iconFA('compress-alt', () => {
      this.showChangingColumns();
      $(showAllColumnsIcon).show();
      $(showChangingColumnsIcon).hide();
    }, 'Hide constant parameters');
    const showAllColumnsIcon = ui.iconFA('expand-arrows-alt', () => {
      this.showAllColumns();
      $(showAllColumnsIcon).hide();
      $(showChangingColumnsIcon).show();
    }, 'Show all parameters');
    $(showChangingColumnsIcon).hide();
    this.setRibbonPanels([...this.getRibbonPanels(), [showChangingColumnsIcon, showAllColumnsIcon]]);

    this.showChangingColumns();

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

    const cache = new DG.LruCache();

    this.grid.onCellPrepare(async (gc) => {
      if (gc.isColHeader || gc.isRowHeader) return;

      if (gc.tableColumn!.name === RUN_NAME_COL_LABEL) {
        gc.customText = '';
        const rows = gc.cell.value.split(' - ') as string[];
        const elems = (rows.length > 1) ? [rows[0], ui.element('br'), rows[1]]: rows;
        gc.element = ui.div(elems, {style: {
          'writing-mode': 'vertical-rl',
          'text-orientation': 'mixed',
          'margin': 'auto',
        }});
        gc.element.parentElement!.style.display = 'flex';
        return;
      }

      if (gc.tableColumn && gc.tableColumn.temp[VIEWER_PATH]) {
        const initialValue = gc.cell.value;

        const viewerConfig = gc.tableColumn.temp[VIEWER_PATH];
        const viewerType = viewerConfig['type'] as string;

        const getElement = async () => {
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
        };

        const uniqueKey = `${initialValue.id}_${viewerType}`;

        if (!cache.get(uniqueKey)) {
          const element = await getElement();
          cache.set(uniqueKey, element);
          gc.element = element;
        } else
          gc.element = cache.get(uniqueKey);

        gc.element.style.width = '100%';
        gc.element.style.height = '100%';
      }
    });
  }
}
