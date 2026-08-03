// ObjectHandler.renderGrid — one place to customize EVERY grid that shows your
// domain rows (the built-in Domain View grid included). The handler decorates a
// grid IN PLACE: column visibility/order, captions, cell renderers.
// Contract: derive presentation from column tags/semTypes only — query frames
// carry no hidden '~item' object column, so never assume one.

const issues = grok.dapi.domains.table('grit.issue');

class IssueGridHandler extends DG.ObjectHandler {
  get type() { return 'grit.issue'; }
  isApplicable(x) { return x instanceof DG.DomainRow && x.typeName === 'grit.issue'; }

  renderGrid(grid, options) {
    grid.columns.byName('description').visible = false;      // hide the long text
    grid.columns.byName('number').move(1);                   // key column first
    (options?.items ?? grid.dataFrame).col('number').setTag(DG.TAGS.FRIENDLY_NAME, '#');
  }
}

// Registration is the standard mechanism — no new API. From now on the Domain
// View grid for grit.issue is decorated by the handler above; a handler WITHOUT
// renderGrid keeps the platform default decoration (tags, ref captions).
DG.ObjectHandler.register(new IssueGridHandler());

// The same handler decorates ad-hoc grids over query results:
const df = await issues.queryDf({limit: 20});
const grid = DG.Grid.create(df);
DG.ObjectHandler.list().find((h) => h.type === 'grit.issue').renderGrid(grid, {items: df});
grok.shell.newView('Issues', [grid.root]);
