/** The portfolio table: the programs of one portfolio (or of every portfolio the user can see) in a
 * Datagrok table view, with the filter panel and the analytics over it: programs by stage and
 * modality, time in stage, the readiness matrix of the next gate and the compounds by role. The
 * matrix and the compound roles are frames of their own, linked to the programs by id, so the
 * filters drive every viewer. Clicking a row opens the program the way the board does. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {COL, PortfolioTables, portfolioTables} from '../portfolio-table';
import {button} from '@datagrok-libraries/u2';
// a cycle with portal-view, safe: neither uses the other at load
import {portalAppView, PortalView} from './portal-view';
import {installStyles} from './styles';
import {openProgramPage} from './my-programs-widget';
import {programAddress} from '../model/pinning';

export const APP_ROOT = '/apps/ProjectPortalCore';
export const APP_PATH = `${APP_ROOT}/Portal`;
/** The last segment of a table view's address. */
export const TABLE_SEGMENT = 'table';

/** `/apps/ProjectPortalCore/Portal/<code>/table`, or `.../Portal/table` for every portfolio. */
export const tablePath = (code?: string | null) => code ? `${APP_PATH}/${code}/${TABLE_SEGMENT}` : `${APP_PATH}/${TABLE_SEGMENT}`;

/** The viewers of the default layout, by title; the BDD scenario and the tests read them. */
export const VIEWER_TITLES = {
  byStage: 'Programs by stage',
  timeInStage: 'Days in stage',
  readiness: 'Readiness for the next gate',
  compounds: 'Compounds by role',
} as const;

const FILTER_COLUMNS = [COL.portfolio, COL.stage, COL.status, COL.modality, COL.ready, COL.readiness,
  COL.daysInStage, COL.attention, COL.workspace, COL.source];

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

/** The current row's program at a glance in the context panel, as a board card shows it; its Open
 * program opens the program's page in the Project Portal. */
function programPanel(tv: DG.TableView, t: PortfolioTables): void {
  installStyles();
  const portal = new PortalView();
  portal.openElsewhere = (code, p) => void openProgramPage({portfolioCode: code, programId: p.id,
    address: programAddress(code, p.code)}).catch((e) => grok.shell.error(e?.message ?? String(e)));
  let loaded: Promise<void> | null = null;
  const df = t.programs;
  const show = async (row: number) => {
    if (row < 0)
      return;
    const id = df.get(COL.id, row) as string;
    const code = df.get(COL.portfolio, row) as string;
    loaded ??= portal.load(code);
    await loaded;
    await portal.previewProgram(code, id);
  };
  tv.subs.push(df.onCurrentRowChanged.subscribe(() => void show(df.currentRowIdx).catch((e) =>
    grok.shell.error(e?.message ?? String(e)))));
}

/** The default layout over [t]: grid, filters, and the four viewers. */
export function layOut(tv: DG.TableView, t: PortfolioTables): void {
  const df = t.programs;
  const several = t.portfolios.length > 1;
  for (const hidden of [COL.id, COL.stageOrdinal, COL.stageInOrder, COL.gaps, ...(several ? [] : [COL.portfolio, COL.setting])])
    tv.grid.columns.byName(hidden)!.visible = false;
  tv.grid.columns.byName(COL.ready)!.cellType = 'bool';
  // a portfolio without attention rules has no attention to show
  if (df.col(COL.attention)!.stats.missingValueCount === df.rowCount) {
    for (const hidden of [COL.attention, COL.attentionRules])
      tv.grid.columns.byName(hidden)!.visible = false;
  }

  // a range filter over a numeric column whose values are all zero drops every row (platform,
  // 2026-10-01): leave such a column out until it varies
  const varies = (c: DG.Column) => !c.isNumerical || c.stats.min < c.stats.max;
  tv.filters({columnNames: FILTER_COLUMNS.filter((c) => df.col(c) && varies(df.col(c)!) && (several || c !== COL.portfolio))});

  const byStage = tv.addViewer(DG.VIEWER.BAR_CHART, {
    title: VIEWER_TITLES.byStage, showTitle: true, splitColumnName: COL.stageInOrder, valueColumnName: COL.program,
    valueAggrType: DG.AGG.TOTAL_COUNT, stackColumnName: COL.modality, barSortType: 'by category', barSortOrder: 'asc',
  });
  tv.dockManager.dock(byStage, DG.DOCK_TYPE.RIGHT, null, VIEWER_TITLES.byStage, 0.33);
  const timeInStage = tv.addViewer(DG.VIEWER.BOX_PLOT, {
    title: VIEWER_TITLES.timeInStage, showTitle: true, categoryColumnNames: [COL.stage], valueColumnName: COL.daysInStage,
  });
  tv.dockManager.dock(timeInStage, DG.DOCK_TYPE.DOWN, tv.dockManager.findNode(byStage.root) ?? null, VIEWER_TITLES.timeInStage, 0.5);

  // the two frames of their own follow the programs' filter through the program id
  grok.data.linkTables(df, t.readiness, [COL.id], [COL.id], [DG.SYNC_TYPE.FILTER_TO_FILTER, DG.SYNC_TYPE.CURRENT_ROW_TO_ROW]);
  grok.data.linkTables(df, t.compounds, [COL.id], [COL.id], [DG.SYNC_TYPE.FILTER_TO_FILTER]);

  // programs × the criteria of their next stage, a criterion's state as its color
  const matrix = DG.Viewer.grid(t.readiness, {title: VIEWER_TITLES.readiness, showTitle: true} as any);
  for (const hidden of [COL.id, ...(several ? [] : [COL.portfolio])])
    matrix.columns.byName(hidden)!.visible = false;
  tv.dockManager.dock(matrix, DG.DOCK_TYPE.DOWN, tv.dockManager.findNode(tv.grid.root) ?? null, VIEWER_TITLES.readiness, 0.4);
  const compounds = DG.Viewer.barChart(t.compounds, {
    title: VIEWER_TITLES.compounds, showTitle: true, splitColumnName: COL.program, valueColumnName: COL.compounds,
    valueAggrType: DG.AGG.SUM, stackColumnName: 'role',
  });
  tv.dockManager.dock(compounds, DG.DOCK_TYPE.DOWN, tv.dockManager.findNode(timeInStage.root) ?? null, VIEWER_TITLES.compounds, 0.4);
  (tv as any).portfolioTables = t;
  (tv as any).extraViewers = [matrix, compounds];
}

/** Opens the table view of [portfolioCode], or of every portfolio the user can see, in the shell. */
export async function portfolioTableView(portfolioCode?: string | null): Promise<DG.TableView> {
  const t = await portfolioTables(portfolioCode ?? undefined);
  const tv = grok.shell.addTableView(t.programs);
  tv.name = portfolioCode ? `${t.portfolios[0].name}: programs` : 'All portfolios: programs';
  // the full address either way: the browser shows it as set. (A view the app function returns
  // also gets the app's root as its basePath, so its `path` getter then reads the root twice.)
  tv.path = tablePath(portfolioCode);
  layOut(tv, t);
  programPanel(tv, t);
  const board = named(button('Board', () =>
    void portalAppView(portfolioCode ?? undefined).then((v) => grok.shell.addView(v))), 'board-button');
  board.title = 'The stage board of this portfolio';
  tv.setRibbonPanels([...tv.getRibbonPanels(), [board]]);
  return tv;
}

/** [portfolioTableView] from a button: a progress indicator, and an error balloon instead of a throw. */
export async function openPortfolioTable(portfolioCode: string | null): Promise<DG.TableView | null> {
  const pi = DG.TaskBarProgressIndicator.create('Reading the portfolio');
  try {
    return await portfolioTableView(portfolioCode);
  } catch (e: any) {
    grok.shell.error(e?.message ?? String(e));
    return null;
  } finally {
    pi.close();
  }
}
