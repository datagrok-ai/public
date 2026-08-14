// DG.DomainView — the platform's Domain View for one table, created
// programmatically: the /domains page's gallery, search, sort, paging, filter
// panel and in-grid editing, either as a view of its own or embedded in yours.

// A full page for the table (what /domains/grit/issue opens).
const issues = DG.DomainView.create({schema: 'grit', table: 'issue'});
grok.shell.addView(issues);
issues.showFilters();                       // server-backed facets and histograms

// The current query — panel facets AND the search string, ready to serialize
// into a deep link or to re-run uncapped through the DomainQuery function.
console.log(issues.query);                  // {schema, table, filters, orderBy, limit}

// Embedded: the same view scoped to a subset, mounted inside your own layout.
const [project] = await grok.dapi.domains.table('grit.project').query({limit: 1});
if (project != null) {
  const scoped = DG.DomainView.create({schema: 'grit', table: 'issue', embedded: true,
    permanentFilter: `project_id = "${project.id}"`});   // invisible, always applied
  grok.shell.newView(`${project.key} issues`, [ui.box(scoped.root)]);
}
