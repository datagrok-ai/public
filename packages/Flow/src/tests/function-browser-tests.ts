/** Regression tests for the function catalog: the Flow-specific exclusion list
 *  and the "what it does" classification. Designed from the live catalog — see
 *  docs/func-catalog-snapshot.md (and the retired func-inventory-diag) for the
 *  data these assertions encode. */
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, EXCLUDED_PACKAGES,
} from '../rete/node-factory';
import {EXCLUDED_FUNC_NQNAMES} from '../rete/excluded-funcs';
import {
  categorizeFunc, FUNC_CATEGORIES, funcMatchesSearch, nameMatchesQuery,
  queryConnectionName, FunctionBrowser, funcOutputsWidget,
} from '../panel/function-browser';
import {getTags} from '../utils/dart-proxy-utils';
import {statusLabel} from '../execution/execution-visualizer';
import {NodeExecStatus} from '../execution/execution-state';

category('Flow: function browser', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('exclusion list keeps the catalog clean', async () => {
    const funcs = getRegisteredFuncs();
    expect(funcs.length > 100, true, 'catalog is non-trivially populated');

    // No dev/test/internal packages.
    const fromExcludedPkg = funcs.filter((f) => EXCLUDED_PACKAGES.includes(f.packageName));
    expect(fromExcludedPkg.length, 0, `excluded packages leaked: ${fromExcludedPkg.map((f) => f.packageName).join(',')}`);

    // No test scaffolding by name.
    const testNamed = funcs.filter((f) => (f.func.name ?? '').toLowerCase().startsWith('test'));
    expect(testNamed.length, 0, `test* funcs leaked: ${testNamed.map((f) => f.func.name).slice(0, 5).join(',')}`);

    // No command/dialog wrappers (funccall inputs).
    const funccallWrappers = funcs.filter((f) => {
      try {return f.func.inputs.some((p) => String(p.propertyType) === 'funccall');} catch {return false;}
    });
    expect(funccallWrappers.length, 0, 'funccall-wrapper funcs leaked');

    // No view-producing functions (a whole TableView can't be composed/previewed).
    const viewOutputs = funcs.filter((f) => {
      try {return f.func.outputs.some((p) => String(p.propertyType) === 'view');} catch {return false;}
    });
    expect(viewOutputs.length, 0, `view-output funcs leaked: ${viewOutputs.map((f) => f.func.name).slice(0, 5).join(',')}`);
  });

  test('denylisted nqNames never survive into the catalog', async () => {
    const survivingNq = new Set(getRegisteredFuncs().map((f) => {
      try {return f.func.nqName;} catch {return f.func.name;}
    }));
    // Every entry that exists on this stand must be filtered out.
    const leaked: string[] = [];
    for (const nq of EXCLUDED_FUNC_NQNAMES)
      if (survivingNq.has(nq)) leaked.push(nq);
    expect(leaked.length, 0, `denylisted funcs leaked: ${leaked.slice(0, 8).join(', ')}`);
  });

  test('machinery tags and right-click actions are excluded (but widgets are kept)', async () => {
    const funcs = getRegisteredFuncs();
    // No surviving func carries a UI-machinery tag — checking tags (not just the
    // role field) is the biggest declutter lever. NOTE: panel/widget/widgets/
    // tooltip are deliberately NOT banned — widget-producing functions are usable
    // in Flow (Widgets pane + preview).
    const banned = new Set(['internal', 'moleculesketcher', 'folderviewer',
      'cellrenderer', 'viewers', 'filehandler', 'semtypedetector', 'apptreebrowser', '@editors']);
    const withBannedTag = funcs.filter((f) => {
      const tokens = getTags(f.func).map((t) => t.trim().toLowerCase());
      return tokens.some((t) => banned.has(t));
    });
    expect(withBannedTag.length, 0,
      `machinery leaked: ${withBannedTag.map((f) => f.func.name).slice(0, 6).join(', ')}`);

    // No semantic_value (right-click) inputs, no filter-DSL-call outputs.
    const semValue = funcs.filter((f) => {
      try {return f.func.inputs.some((p) => String(p.propertyType) === 'semantic_value');} catch {return false;}
    });
    expect(semValue.length, 0, 'semantic_value right-click actions leaked');
    const filterCalls = funcs.filter((f) => {
      try {
        return f.func.outputs.some((p) =>
          ['tablerowfiltercall', 'colfiltercall'].includes(String(p.propertyType)));
      } catch {return false;}
    });
    expect(filterCalls.length, 0, 'filter-DSL builder funcs leaked');
  });

  test('widget-producing functions are kept (Widgets pane populated)', async () => {
    // Widgets are supported (preview) — a function that outputs a widget must not
    // be excluded just for being a widget/panel. There should be widget nodes.
    const widgets = getRegisteredFuncs().filter(funcOutputsWidget);
    expect(widgets.length > 0, true, 'at least one widget-producing function survives');
    // And a right-click widget (semantic_value input) is still excluded.
    const ctxWidgets = widgets.filter((f) => {
      try {return f.func.inputs.some((p) => String(p.propertyType) === 'semantic_value');} catch {return false;}
    });
    expect(ctxWidgets.length, 0, 'context (semantic_value) widgets are still excluded');
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

  test('the toolbox renders a Cheminformatics section when chem funcs exist', async () => {
    if (!getRegisteredFuncs().some((f) => f.packageName === 'Chem')) {
      expect(true, true, 'no Chem funcs on this stand — skipped');
      return;
    }
    const browser = new FunctionBrowser({
      onFunctionDoubleClick: () => {}, onBuiltinNodeDoubleClick: () => {}, onFileDoubleClick: () => {},
    });
    document.body.appendChild(browser.root);
    try {
      browser.render();
      const header = browser.root.querySelector('[data-testid="ff-browser-section-cheminformatics"]');
      expect(!!header, true, 'Cheminformatics section header present');
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
      onFunctionDoubleClick: () => {}, onBuiltinNodeDoubleClick: () => {}, onFileDoubleClick: () => {},
    });
    document.body.appendChild(browser.root);
    try {
      browser.render();

      // The Files pane (open by default) and the Queries pane both exist.
      expect(!!browser.root.querySelector('[data-testid="ff-browser-files"]'), true, 'Files pane present');
      const queriesPane = browser.root.querySelector('[data-testid="ff-browser-queries"]') as HTMLElement | null;
      expect(!!queriesPane, true, 'Queries pane present');

      // At least one per-connection sub-section, each tagged with its connection.
      const connSections = browser.root.querySelectorAll('[data-testid^="ff-browser-query-conn"]');
      expect(connSections.length > 0, true, 'queries split into per-connection sub-sections');
      expect((connSections[0] as HTMLElement).dataset.queryConn != null, true, 'sub-section carries data-query-conn');

      // A known query lives INSIDE the Queries pane and NOT in any category section.
      const sample = queries[0];
      const items = Array.from(browser.root.querySelectorAll(`[data-func="${sample.func.name}"]`)) as HTMLElement[];
      expect(items.length > 0, true, `query ${sample.func.name} appears in the toolbox`);
      for (const it of items)
        expect(queriesPane!.contains(it), true, `query item ${sample.func.name} is under the Queries pane`);
    } finally {
      browser.root.remove();
    }
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
