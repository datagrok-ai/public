/** Regression tests for the function catalog: the Flow-specific exclusion list
 *  and the "what it does" classification. Designed from the live catalog — see
 *  docs/func-catalog-snapshot.md (and the retired func-inventory-diag) for the
 *  data these assertions encode. */
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, EXCLUDED_PACKAGES,
} from '../rete/node-factory';
import {categorizeFunc, FUNC_CATEGORIES} from '../panel/function-browser';
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
      const cat = categorizeFunc(info.func, info.role);
      expect((FUNC_CATEGORIES as readonly string[]).includes(cat), true, `unknown category ${cat}`);
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
