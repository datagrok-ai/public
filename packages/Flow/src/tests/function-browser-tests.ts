/** Regression tests for the function catalog: the allowlist-based inclusion
 *  rule and the "what it does" classification. Designed from the live catalog —
 *  see docs/func-catalog-snapshot.md for the data these assertions encode. */
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, isWorkflowFunc, shouldIncludeFunc,
} from '../rete/node-factory';
import {INCLUDED_FUNC_NQNAMES} from '../rete/included-funcs';
import {
  categorizeFunc, FUNC_CATEGORIES, funcMatchesSearch, nameMatchesQuery,
  queryConnectionName, FunctionBrowser, funcOutputsWidget, orderDomainSection,
} from '../panel/function-browser';
import type {FuncInfo} from '../rete/node-factory';
import {statusLabel} from '../execution/execution-visualizer';
import {NodeExecStatus} from '../execution/execution-state';

category('Flow: function browser', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('the catalog is allowlist-based: only explicitly included funcs survive', async () => {
    const funcs = getRegisteredFuncs();
    expect(funcs.length > 100, true, 'catalog is non-trivially populated');

    // Known allowlisted core funcs are present.
    const names = new Set(funcs.map((f) => f.func.name));
    for (const known of ['JoinTables', 'OpenFile', 'AddNewColumn', 'Aggregate'])
      expect(names.has(known), true, `allowlisted ${known} present`);

    // Formerly-denied machinery is NOT on the allowlist and stays out: dev/test
    // packages, denylist-era helpers, panels/sketchers — none of them listed.
    for (const gone of ['Chem:getRdKitModule', 'core:BatchCall', 'core:Project'])
      expect(INCLUDED_FUNC_NQNAMES.has(gone), false, `${gone} not on the include list`);
    const devPkgs = new Set(['Dbtests', 'ApiTests', 'UiTests', 'DevTools', 'Tutorials', 'ApiSamples', 'UsageAnalysis']);
    const fromDevPkg = funcs.filter((f) => devPkgs.has(f.packageName));
    expect(fromDevPkg.length, 0, `dev/test packages leaked: ${fromDevPkg.map((f) => f.packageName).join(',')}`);
  });

  test('shouldIncludeFunc: allowlist, opt-in, opt-out precedence', async () => {
    const fake = (nqName: string, options: Record<string, unknown> = {}): DG.Func =>
      ({nqName, name: nqName.split(':').pop(), options} as unknown as DG.Func);
    expect(shouldIncludeFunc(fake('core:JoinTables')), true, 'allowlisted → in');
    expect(shouldIncludeFunc(fake('SomePkg:notListedAnywhere')), false, 'unlisted → out');
    expect(shouldIncludeFunc(fake('SomePkg:optedIn', {includeInFlow: 'true'})), true, 'meta opt-in (string) → in');
    expect(shouldIncludeFunc(fake('SomePkg:optedIn2', {includeInFlow: true})), true, 'meta opt-in (bool) → in');
    expect(shouldIncludeFunc(fake('core:JoinTables', {includeInFlow: 'false'})), false,
      'meta opt-out beats the allowlist');
    // A saved flow (Script with language `flow`) is always included, listed or not.
    const flowScript = DG.Script.create('//name: MyFlow\n//language: flow\n');
    expect(shouldIncludeFunc(flowScript), true, 'workflow → always in');
    // A query (DG.DataQuery) is always included by kind — never needs listing.
    const query = getRegisteredFuncs().map((f) => f.func).find((f) => f instanceof DG.DataQuery);
    if (query) {
      expect(shouldIncludeFunc(query), true, 'query → always in');
      expect(INCLUDED_FUNC_NQNAMES.has(query.nqName), false, 'queries are not on the static list');
    }
  });

  test('meta.includeInFlow: false opts a function out of the catalog', async () => {
    // No surviving function declares the opt-out (meta.* surfaces as
    // func.options; the value may arrive as boolean false or the string 'false').
    const leaked = getRegisteredFuncs().filter((f) => {
      try {
        const v = (f.func.options as Record<string, unknown>)?.['includeInFlow'];
        return v === false || String(v).toLowerCase() === 'false';
      } catch {return false;}
    });
    expect(leaked.length, 0, `includeInFlow:false funcs leaked: ${leaked.map((f) => f.func.name).join(', ')}`);

    // End-to-end: Flow's own dialog opener declares the opt-out in its decorator
    // meta — it must exist on the stand yet be absent from the catalog.
    const exists = DG.Func.find({name: 'openCreationScriptFlowDialog'}).length > 0;
    if (!exists) return; // older Flow build on this stand — skip the e2e half
    const inCatalog = getRegisteredFuncs().some((f) => f.func.name === 'openCreationScriptFlowDialog');
    expect(inCatalog, false, 'openCreationScriptFlowDialog is opted out via meta.includeInFlow');
  });

  test('widget-producing functions are kept (Widgets pane populated)', async () => {
    // Widgets are supported (preview) — allowlisted widget-producing functions
    // survive into the Widgets pane.
    const widgets = getRegisteredFuncs().filter(funcOutputsWidget);
    expect(widgets.length > 0, true, 'at least one widget-producing function survives');
    // And right-click (semantic_value) widgets were never allowlisted.
    const ctxWidgets = widgets.filter((f) => {
      try {return f.func.inputs.some((p) => String(p.propertyType) === 'semantic_value');} catch {return false;}
    });
    expect(ctxWidgets.length, 0, 'context (semantic_value) widgets stay out');
  });

  test('categorizeFunc places funcs by what they do', async () => {
    // (functionName -> expected category). Guarded: only assert when the func
    // exists on this stand, so the test is portable across deployments.
    const cases: Array<[string, string]> = [
      ['JoinTables', 'Combine Tables'],   // 2 table inputs -> NOT a data source
      ['LinkTables', 'Combine Tables'],
      ['OpenFile', 'Data Sources'],       // table out, no table in
      ['OpenTable', 'Data Sources'],
      ['Aggregate', 'Transform Tables'],  // 1 table in -> table out
      ['Unpivot', 'Transform Tables'],
      ['AddNewColumn', 'Column Operations'],
    ];
    let checked = 0;
    for (const [name, expected] of cases) {
      const f = DG.Func.find({name})[0];
      if (!f) continue;
      checked++;
      expect(categorizeFunc(f, null), expected, `${name} -> ${expected}`);
    }
    expect(checked > 0, true, 'at least one known function was available to classify');
  });

  test('saved flows classify as Workflows in every grouping', async () => {
    // A flow entity is a DG.Script with language `flow`; its signature would
    // otherwise read as a data source — it gets its own section instead.
    const flowScript = DG.Script.create('//name: MyFlow\n//language: flow\n');
    expect(isWorkflowFunc(flowScript), true, 'a flow script is a workflow');
    expect(categorizeFunc(flowScript, null, 'Flow'), 'Workflows', 'category routing');
    expect((FUNC_CATEGORIES as readonly string[]).includes('Workflows'), true,
      'Workflows is an ordered category');

    const jsScript = DG.Script.create('//name: MyJs\n//language: javascript\nlet x = 1;');
    expect(isWorkflowFunc(jsScript), false, 'other script languages stay in their signature category');
    const plainFunc = DG.Func.find({name: 'OpenFile'})[0];
    if (plainFunc) expect(isWorkflowFunc(plainFunc), false, 'a non-script func is not a workflow');
  });

  test('Data Sources leads the category order', async () => {
    expect(FUNC_CATEGORIES[0], 'Data Sources');
    // Every classified func lands in a known category.
    for (const info of getRegisteredFuncs().slice(0, 200)) {
      const cat = categorizeFunc(info.func, info.role, info.packageName);
      expect((FUNC_CATEGORIES as readonly string[]).includes(cat), true, `unknown category ${cat}`);
    }
  });

  test('chem/bio operations group into their domain sections', async () => {
    const hasDataInput = (f: {func: DG.Func}): boolean => {
      try {return f.func.inputs.some((p) => ['dataframe', 'column', 'column_list'].includes(String(p.propertyType)));}
      catch {return false;}
    };
    // A Chem/Bio function that OPERATES on data (dataframe/column input) lands in
    // Cheminformatics/Bioinformatics — the domain wins over the task category.
    const chem = getRegisteredFuncs().find((f) => f.packageName === 'Chem' && hasDataInput(f));
    if (chem) expect(categorizeFunc(chem.func, chem.role, 'Chem'), 'Cheminformatics');
    const bio = getRegisteredFuncs().find((f) => f.packageName === 'Bio' && hasDataInput(f));
    if (bio) expect(categorizeFunc(bio.func, bio.role, 'Bio'), 'Bioinformatics');
    // A core/general func is unaffected by the domain routing.
    const join = DG.Func.find({name: 'JoinTables'})[0];
    if (join) expect(categorizeFunc(join, null, ''), 'Combine Tables');

    // The two domain sections are in the ordered category list, after Visualize.
    expect((FUNC_CATEGORIES as readonly string[]).includes('Cheminformatics'), true);
    expect((FUNC_CATEGORIES as readonly string[]).includes('Bioinformatics'), true);
    expect(FUNC_CATEGORIES.indexOf('Cheminformatics') > FUNC_CATEGORIES.indexOf('Visualize'), true);
  });

  test('chem/bio *sources* (no data input) stay out of the domain sections', async () => {
    // A pure source/query/generator (produces a table from scalars, no dataframe/
    // column input) is NOT an operation — it falls back to its task category, so
    // the domain sections hold only functions that do something to your data.
    const chemSource = getRegisteredFuncs().find((f) => {
      if (f.packageName !== 'Chem' && f.packageName !== 'Chembl' && f.packageName !== 'ChemblApi') return false;
      try {
        const ins = f.func.inputs.map((p) => String(p.propertyType));
        const outs = f.func.outputs.map((p) => String(p.propertyType));
        return outs.includes('dataframe') && !ins.some((t) => ['dataframe', 'column', 'column_list'].includes(t));
      } catch {return false;}
    });
    if (!chemSource) {
      expect(true, true, 'no chem source on this stand — skipped');
      return;
    }
    const cat = categorizeFunc(chemSource.func, chemSource.role, chemSource.packageName);
    expect(cat !== 'Cheminformatics', true, `${chemSource.func.name} is a source, not a chem operation (got ${cat})`);
  });

  test('the toolbox floats Cheminformatics/Bioinformatics to the top of the categories', async () => {
    if (!getRegisteredFuncs().some((f) => f.packageName === 'Chem')) {
      expect(true, true, 'no Chem funcs on this stand — skipped');
      return;
    }
    const browser = new FunctionBrowser({
      onFunctionDoubleClick: () => {}, onBuiltinNodeDoubleClick: () => {},
      onFileDoubleClick: () => {}, onLocalFilesPicked: () => {},
    });
    document.body.appendChild(browser.root);
    try {
      browser.render();
      const chem = browser.root.querySelector('[data-testid="ff-browser-section-cheminformatics"]');
      expect(!!chem, true, 'Cheminformatics section header present');
      // Order (Files/Queries/Workflows live in the top tabs now): domain
      // sections → task categories (Data Sources…) → Viewers → built-ins
      // (Inputs…) → Other → Debug last.
      const viewers = browser.root.querySelector('[data-testid="ff-browser-viewers"]');
      const inputs = browser.root.querySelector('[data-testid="ff-browser-section-inputs"]');
      const dataSources = browser.root.querySelector('[data-testid="ff-browser-section-data-sources"]');
      const other = browser.root.querySelector('[data-testid="ff-browser-section-other"]');
      const debug = browser.root.querySelector('[data-testid="ff-browser-section-debug"]');
      const before = (a: Element | null, b: Element | null): boolean =>
        !!a && !!b && !!(a.compareDocumentPosition(b) & Node.DOCUMENT_POSITION_FOLLOWING);
      if (dataSources) expect(before(chem, dataSources), true, 'Cheminformatics comes before Data Sources');
      if (viewers && dataSources) expect(before(dataSources, viewers), true, 'task categories come before Viewers');
      if (viewers && inputs) expect(before(viewers, inputs), true, 'Viewers comes before the Inputs built-ins');
      // Debug is a static built-in (Breakpoint) — it exists on every stand.
      expect(!!debug, true, 'Debug section present');
      if (other) expect(before(inputs, other), true, 'Other comes after the built-ins');
      if (other) expect(before(other, debug), true, 'Debug comes last, after Other');
      // The accordion holds no Files/Queries/Workflows panes anymore.
      expect(browser.accordion!.panes.some((p) => ['Files', 'Queries', 'Workflows'].includes(p.name)),
        false, 'collection panes moved out of the accordion');
    } finally {
      browser.root.remove();
    }
  });

  test('search matches the raw func name even when the toolbox shows friendlyName', async () => {
    // The bug: "OpenFile" returned nothing because the list shows "Open File".
    const info = getRegisteredFuncs().find((f) => f.func.name === 'OpenFile');
    if (!info) {
      expect(true, true, 'OpenFile not on this stand — skipped');
      return;
    }
    expect(info.name !== info.func.name, true, 'display name differs from func name (Open File vs OpenFile)');
    expect(funcMatchesSearch(info, 'openfile'), true, 'matches the raw function name');
    expect(funcMatchesSearch(info, 'open file'), true, 'still matches the friendly/display name');
    expect(funcMatchesSearch(info, 'open'), true, 'matches a shared prefix');
    expect(funcMatchesSearch(info, 'zzzznope'), false, 'no false positive');
  });

  test('nameMatchesQuery is case- and whitespace-insensitive', async () => {
    // Built-in items (e.g. "Table Output") are matched by display name, so the
    // matcher must accept both spaced and unspaced queries.
    expect(nameMatchesQuery('Table Output', 'tableoutput'), true, 'unspaced query');
    expect(nameMatchesQuery('Table Output', 'table output'), true, 'spaced query');
    expect(nameMatchesQuery('Open File', 'openfile'), true, 'unspaced');
    expect(nameMatchesQuery('Open File', 'OPEN'), true, 'case-insensitive prefix');
    expect(nameMatchesQuery('Table Output', 'value'), false, 'no false positive');
  });

  test('queries are grouped by connection and kept out of the categories', async () => {
    const queries = getRegisteredFuncs().filter((f) => f.func instanceof DG.DataQuery);
    if (queries.length === 0) {
      expect(true, true, 'no queries on this stand — skipped');
      return;
    }
    // Every query reports a non-empty connection name (the grouping key).
    for (const q of queries.slice(0, 20))
      expect(queryConnectionName(q).length > 0, true, `query ${q.func.name} has a connection name`);

    const browser = new FunctionBrowser({
      onFunctionDoubleClick: () => {}, onBuiltinNodeDoubleClick: () => {},
      onFileDoubleClick: () => {}, onLocalFilesPicked: () => {},
    });
    document.body.appendChild(browser.root);
    try {
      browser.render();

      // Queries live in the top Queries TAB now — activate it so its content
      // attaches, then expand every per-connection sub-pane so items
      // materialize in the DOM.
      browser.showTab('Queries');
      const queriesPane = browser.root.querySelector('[data-testid="ff-browser-queries"]') as HTMLElement | null;
      expect(!!queriesPane, true, 'Queries tab content present');
      const connSections = browser.root.querySelectorAll('[data-testid^="ff-browser-query-conn"]');
      expect(connSections.length > 0, true, 'queries split into per-connection sub-sections');
      expect((connSections[0] as HTMLElement).dataset.queryConn != null, true, 'sub-section carries data-query-conn');
      for (const p of browser.queriesAccordion!.panes) p.expanded = true;

      // A known query lives INSIDE the Queries tab and NOT in any category section.
      const sample = queries[0];
      const items = Array.from(browser.root.querySelectorAll(`[data-func="${sample.func.name}"]`)) as HTMLElement[];
      expect(items.length > 0, true, `query ${sample.func.name} appears in the toolbox`);
      for (const it of items)
        expect(queriesPane!.contains(it), true, `query item ${sample.func.name} is under the Queries tab`);
    } finally {
      browser.root.remove();
    }
  });

  test('toolbox sections are platform accordion panes with self-persisted state', async () => {
    const browser = new FunctionBrowser({
      onFunctionDoubleClick: () => {}, onBuiltinNodeDoubleClick: () => {},
      onFileDoubleClick: () => {}, onLocalFilesPicked: () => {},
    });
    document.body.appendChild(browser.root);
    const lsKey = 'Accordion:funcflow.toolbox';
    const saved = localStorage.getItem(lsKey);
    try {
      browser.render();
      // The sections ARE a DG.Accordion — no custom collapsible divs left.
      expect(!!browser.root.querySelector('.d4-accordion'), true, 'platform accordion present');
      expect(browser.root.querySelector('.funcflow-section-header') == null, true, 'no custom section headers');
      // Guide/test hooks survive on the pane headers.
      const inputsHeader = browser.root.querySelector(
        '[data-testid="ff-browser-section-inputs"]') as HTMLElement | null;
      expect(!!inputsHeader, true, 'Inputs pane header carries its test id');
      expect(inputsHeader!.dataset.section, 'Inputs', 'header keeps data-section for the guide');

      // Clicking a header persists the pane state under the accordion key
      // (the platform writes localStorage["Accordion:<key>"] itself)...
      expect(browser.accordion!.getPane('Inputs').expanded, false, 'Inputs starts collapsed');
      inputsHeader!.click();
      expect(browser.accordion!.getPane('Inputs').expanded, true, 'click expands the pane');
      const stored = JSON.parse(localStorage.getItem(lsKey) ?? '{}') as Record<string, boolean>;
      expect(stored['Inputs'], true, 'expanded state persisted by the platform');

      // ...and a fresh render restores it.
      browser.render();
      expect(browser.accordion!.getPane('Inputs').expanded, true, 'state restored on re-render');
    } finally {
      browser.root.remove();
      if (saved == null) localStorage.removeItem(lsKey);
      else localStorage.setItem(lsKey, saved);
    }
  });

  test('orderDomainSection floats the flagship package + most-used ops to the top', async () => {
    // Only name / func.name / packageName drive the sort.
    const mk = (name: string, pkg: string): FuncInfo => ({
      func: {name} as DG.Func, name, role: null, tags: [], packageName: pkg,
      nodeTypeName: `DG Functions/Cheminformatics/${name}`,
    });
    // Deliberately shuffled: a non-Chem pkg, a plain Chem func, and Chem's
    // descriptors/properties which must lead.
    const items = [
      mk('Chemspace Search', 'Chemspace'),
      mk('Add Chem Risks', 'Chem'),
      mk('addChemPropertiesColumns', 'Chem'),
      mk('getMorganFingerprints', 'Chem'),
      mk('Admetica Eval', 'Admetica'),
      mk('Get Descriptors', 'Chem'),
    ];
    const ordered = orderDomainSection(items, 'Cheminformatics').map((i) => i.packageName);
    // Every Chem item precedes every non-Chem item.
    const lastChem = ordered.lastIndexOf('Chem');
    const firstOther = ordered.findIndex((p) => p !== 'Chem');
    expect(lastChem < firstOther, true, 'all Chem funcs come before other cheminformatics pkgs');

    // Within Chem, descriptors/properties/fingerprints lead the plain "Add Chem Risks".
    const names = orderDomainSection(items, 'Cheminformatics').map((i) => i.name);
    const idx = (n: string): number => names.indexOf(n);
    expect(idx('Get Descriptors') < idx('Add Chem Risks'), true, 'descriptors before a generic risks func');
    expect(idx('addChemPropertiesColumns') < idx('Add Chem Risks'), true, 'properties before generic');
    expect(idx('getMorganFingerprints') < idx('Add Chem Risks'), true, 'fingerprints before generic');

    // Bio flagship routing works too.
    const bioItems = [mk('Helm Convert', 'Helm'), mk('Sequence Descriptors', 'Bio')];
    const bioOrdered = orderDomainSection(bioItems, 'Bioinformatics').map((i) => i.packageName);
    expect(bioOrdered[0], 'Bio', 'Bio leads the Bioinformatics section');
  });

  test('statusLabel renders plain-language status', async () => {
    expect(statusLabel(NodeExecStatus.idle), '');
    expect(statusLabel(NodeExecStatus.running), 'Running…');
    expect(statusLabel(NodeExecStatus.completed), 'Done');
    expect(statusLabel(NodeExecStatus.completed, '1,204 × 8'), 'Done · 1,204 × 8');
    expect(statusLabel(NodeExecStatus.errored), 'Error');
    expect(statusLabel(NodeExecStatus.stale), 'Out of date');
  });
});
