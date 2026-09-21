import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, before, test, after, expect, delay} from '@datagrok-libraries/test/src/test';
import {TestManager} from '../package-testing';


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
    await delay(100);
    expect(testManager.testsResultsDf.get('package', 0), 'DevTools');
    expect(testManager.testsResultsDf.get('category', 0), 'FSE');
    expect(testManager.testsResultsDf.get('name', 0), 'exist');
  });


  after(async () => {
    testManager.testManagerView.close();
  });
});
