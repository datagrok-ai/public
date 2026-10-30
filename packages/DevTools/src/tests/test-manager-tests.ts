import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, before, test, after, expect, waitUntil} from '@datagrok-libraries/test/src/test';
import {NODE_TYPE, TestManager} from '../package-testing';


category('Test manager', () => {
  const viewName = 'test_manager_test';
  let testManager: TestManager;

  before(async () => {
    testManager = new TestManager(viewName);
    await testManager.init();
  });

  test('Test manager opens', async () => {
    let viewCreated = false;
    for (const v of grok.shell.views) {
      if (v.name === viewName)
        viewCreated = true;
    }
    expect(viewCreated, true);
    expect(testManager.testFunctions.length > 0, true);
  });

  test('Tests are running', async () => {
    const devToolsNode = testManager.tree.items.filter((it) => it.text === 'Dev Tools')[0] as DG.TreeViewGroup;
    const f = testManager.testFunctions.filter((it) => it.package.name === 'DevTools')[0];
    await testManager.collectPackageTests(devToolsNode, f);
    // TreeViewNode.text is the node's text, not its rendered HTML — matching markup here
    // silently selected nothing, and runTestsForSelectedNode() then returned without running.
    testManager.selectedNode = testManager.tree.items.filter((it) => it.text === 'exist')[0];
    expect(testManager.selectedNode != null, true);
    await testManager.runTestsForSelectedNode();
    await waitUntil(() => testManager.testsResultsDf?.rowCount > 0, 'the run left no result');
    expect(testManager.testsResultsDf.get('package', 0), 'DevTools');
    expect(testManager.testsResultsDf.get('category', 0), 'FSE');
    expect(testManager.testsResultsDf.get('name', 0), 'exist');
  });

  test('Unhandled exception row', async () => {
    const condition = {package: 'DevTools', category: 'FSE', name: 'exist'};
    const testInfo = testManager.getTestsInfoGrid(condition, NODE_TYPE.TEST, false, 'unhandled error').testInfo;
    const last = testInfo.rowCount - 1;
    expect(testInfo.get('category', last), 'unhandled');
    expect(testInfo.get('success', last), false);
    expect(testInfo.get('result', last), 'unhandled error');
    expect(testInfo.get('package', last), 'DevTools');
  });

  // Older test libraries return no 'package' column (and, before 2024, no 'logs')
  test('Tests details after running old test library tests', async () => {
    const call = grok.functions.call;
    const savedResults = testManager.testsResultsDf;
    testManager.testsResultsDf = undefined;
    grok.functions.call = async (name: string, params?: {}) => {
      const res = await call.call(grok.functions, name, params);
      if (name === 'DevTools:test') {
        for (const col of ['logs', 'package'])
          if (res.col(col))
            res.columns.remove(col);
      }
      return res;
    };
    try {
      testManager.selectedNode = testManager.tree.items.filter((it) => it.text === 'exist')[0];
      await testManager.runTestsForSelectedNode();
      const condition = {package: 'DevTools', category: 'FSE', name: 'exist'};
      const testInfo = testManager.getTestsInfoGrid(condition, NODE_TYPE.TEST).testInfo;
      expect(testInfo.rowCount, 1);
    }
    finally {
      grok.functions.call = call;
      testManager.testsResultsDf = savedResults;
    }
  });


  after(async () => {
    testManager.testManagerView.close();
  });
});
