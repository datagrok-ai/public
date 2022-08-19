import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Test } from '@datagrok-libraries/utils/src/test';
import { c, _package } from './package';
import { Menu } from 'datagrok-api/dg';

enum NODE_TYPE {
  PACKAGE = 'Package',
  CATEGORY = 'Category',
  TEST = 'TEST'
}

interface ITestManagerUI {
  testsTree: DG.TreeViewNode;
  runButton: HTMLButtonElement;
  runAllButton: HTMLButtonElement;
}

interface IPackageTests {
  name: string;
  friendlyName: string;
  categories: {[index: string]: ICategory};
  totalTests: number;
}

interface ICategory {
  name: string,
  fullName: string,
  tests: IPackageTest[],
  subcategories: {[index: string]: ICategory},
  packageName: string,
  subCatsNames: string[],
  resultDiv: HTMLElement,
  expand: boolean,
  totalTests: number,
}

interface IPackageTest {
  test: Test,
  packageName: string,
  resultDiv: HTMLElement,
}

interface ITestFromUrl {
  packName: string,
  catName: string,
  testName: string
}


let packagesTests: IPackageTests[] = [];
let testsResultsDf: DG.DataFrame;
let testFunctions: any[];
let testManagerView: DG.View;
let selectedNode: DG.TreeViewGroup|DG.TreeViewNode;
const nodeDict: {[id: string]: any} = {};
let debugMode = false;
let benchmarkMode = false;

export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  grok.shell.addView(view);
  return view;
}

export async function createTestManagerView(): Promise<void> {
  let pathSegments = window.location.pathname.split('/');
  pathSegments = pathSegments.map(it => it ? it.replace(/%20/g, ' ') : undefined);
  testFunctions = await collectPackages();
  testManagerView = DG.View.create();
  const testFromUrl = pathSegments.length > 4 ? { packName: pathSegments[4], catName: pathSegments[5], testName: pathSegments[6]} : null;
  const testUIElements: ITestManagerUI = await createTestManagerUI(testFromUrl);
  testManagerView.name = 'Test Manager'
  addView(testManagerView);
  testManagerView.temp['ignoreCloseAll'] = true;
  testManagerView.setRibbonPanels(
    [
      [testUIElements.runButton], 
      [testUIElements.runAllButton], 
      [ui.boolInput('Debug', false, () => { debugMode = !debugMode; }).root],
      [ui.switchInput('Benchmark', false, () => { benchmarkMode = !benchmarkMode; }).root]
    ],
  );
  testManagerView.append(testUIElements.testsTree.root);
  runTestsForSelectedNode();
}

async function collectPackages(packageName?: string): Promise<any[]>  {
  let testFunctions = DG.Func.find({ name: 'Test' , meta: {file: 'package-test.js'}});
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  return testFunctions;
}

async function collectPackageTests(packageNode: DG.TreeViewGroup, f: any, testFromUrl?: ITestFromUrl) {
  if (testFunctions.filter(it => it.package.name === f.package.name).length !== 0 &&
    packagesTests.filter(pt => pt.name === f.package.name).length === 0) {
    await f.package.load({ file: f.options.file });
    const testModule = f.package.getModule(f.options.file);
    if (!testModule) {
      console.error(`Error getting tests from '${f.package.name}/${f.options.file}' module.`);
    }
    const allPackageTests = testModule ? testModule.tests : undefined;
    let packageTestsFinal: { [cat: string]: ICategory } = {};
    if (allPackageTests) {
      Object.keys(allPackageTests).forEach((cat) => {
        const fullCatName = cat;
        const tests: IPackageTest[] = allPackageTests[cat].tests.map((t) => {
          return { test: t, active: true, packageName: f.package.name, fullCatName: fullCatName };
        });
        const subcats = cat.split(':');
        let subcatsFromUrl = [];
        if (testFromUrl && testFromUrl.catName) {
          subcatsFromUrl = testFromUrl.catName.split(':');
        }
        let previousCat = packageTestsFinal;
        for (let i = 0; i < subcats.length; i++) {
          if (!packageTestsFinal[subcats[i]]) {
            previousCat[subcats[i]] = { 
              tests: [], 
              subcategories: {}, 
              packageName: f.package.name,
              name: subcats[i],
              fullName: subcats.slice(0, i+1).join(':'), 
              subCatsNames: [cat],
              resultDiv: null,
              expand: subcats[i] === subcatsFromUrl[i],
              totalTests: 0};
          } else {
            previousCat[subcats[i]].subCatsNames.push(cat);
          };
          if (i === subcats.length - 1) {
            previousCat[subcats[i]].tests = previousCat[subcats[i]].tests ? previousCat[subcats[i]].tests.concat(tests) : tests;
          } else {
            previousCat = previousCat[subcats[i]].subcategories;
          }
        }
      });
      let testsNumInPack = 0;
      Object.keys(packageTestsFinal).forEach(cat => {
        testsNumInPack += addCategoryRecursive(packageNode, packageTestsFinal[cat], testFromUrl);
      })
      packagesTests.push({
        name: f.package.name,
        friendlyName: f.package.friendlyName,
        categories: packageTestsFinal,
        totalTests: testsNumInPack,
      })
    };
  }
}

function addCategoryRecursive(node: DG.TreeViewGroup, category: ICategory, testFromUrl?: ITestFromUrl) {
  const testPassed = ui.div();
  category.resultDiv = testPassed;
  const subnode = node.group(category.name, null, category.expand);
  subnode.root.children[0].append(testPassed);
  setRunTestsMenuAndLabelClick(subnode, category, NODE_TYPE.CATEGORY);
  if (testFromUrl && testFromUrl.catName === category.fullName) {
    selectedNode = subnode;
  }
  const subcats = Object.keys(category.subcategories);
  if (subcats.length > 0) {
    subcats.forEach((subcat) => {
      category.totalTests += addCategoryRecursive(subnode, category.subcategories[subcat], testFromUrl);
    })
  }
  category.tests.forEach(t => {
    const testPassed = ui.div();
    const itemDiv = ui.divH([
      testPassed,
      ui.divText(t.test.name),
    ]);
    const item = subnode.item(itemDiv);
    t.resultDiv = testPassed;
    setRunTestsMenuAndLabelClick(item, t, NODE_TYPE.TEST);
    ui.tooltip.bind(item.root, () => getTestsInfoGrid( `Package = ${t.packageName} and category = ${t.test.category} and name =  ${t.test.name}`, NODE_TYPE.TEST));
    if (testFromUrl && testFromUrl.catName === category.fullName && testFromUrl.testName === t.test.name) {
      selectedNode = item;
    }
  });
  category.totalTests += category.tests.length;
  return category.totalTests
}


async function createTestManagerUI(testFromUrl: ITestFromUrl): Promise<ITestManagerUI> {
  const tree = ui.tree();
  tree.onSelectedNodeChanged.subscribe((res) => {
    selectedNode = res;
  });
  tree.root.style.width = '100%';
  tree.root.addEventListener('keyup', async (e) => {
    if(e.key === 'Enter') {
      runTestsForSelectedNode();
    }
  });
  for (let pack of testFunctions) {
    const testPassed = ui.div();
    pack.resultDiv = testPassed;
    const packNode = tree.group(pack.package.friendlyName, null, testFromUrl && pack.package.name === testFromUrl.packName);
    setRunTestsMenuAndLabelClick(packNode, pack, NODE_TYPE.PACKAGE);
    if (testFromUrl && testFromUrl.packName === pack.package.name) {
      selectedNode = packNode;
      await collectPackageTests(packNode, pack, testFromUrl);
    }
    packNode.root.children[0].append(testPassed);
    packNode.onNodeExpanding.subscribe(() => {
      collectPackageTests(packNode, pack);
    });
  }

  const runTestsButton = ui.bigButton('Run', async () => {
    runTestsForSelectedNode();
  });

  const runAllButton = ui.bigButton('Run All', async () => {
    const nodes = tree.items;
    for (let node of nodes) {
      selectedNode = node;
      await runTestsForSelectedNode();
    }
  });

  return {runButton: runTestsButton, runAllButton: runAllButton, testsTree: tree};
}

async function runTestsForSelectedNode() {
  if (selectedNode) {
    const id = selectedNode.root.id;
      await runAllTests(selectedNode, nodeDict[id].tests, nodeDict[id].nodeType);
    }
}

function setRunTestsMenuAndLabelClick(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
  node.captionLabel.addEventListener('contextmenu', (e) => {
    Menu.popup()
      .item('Run', async () => {
        runAllTests(node, tests, nodeType);
      })
      .show();
    e.preventDefault();
    e.stopPropagation();
  });
  node.captionLabel.addEventListener('click', () => {
    grok.shell.o = getTestsInfoPanel(node, tests, nodeType);
  });
  if (nodeType === NODE_TYPE.TEST) {
    node.captionLabel.addEventListener('dblclick', async () => {
          runAllTests(node, tests, nodeType);
    });
  }
    const nodeId = Object.keys(nodeDict).length.toString();
    node.root.id = nodeId;
    nodeDict[Object.keys(nodeDict).length] = {tests: tests, nodeType: nodeType};
}


function updateTestResultsIcon (resultDiv: HTMLElement, success?: boolean) {
  success === undefined ? resultDiv.innerHTML = '' : updateIcon(success, resultDiv);
}

function testInProgress (resultDiv: HTMLElement, running: boolean) {
  let icon = ui.iconFA('spinner-third');
  icon.classList.add('fa-spin');
  icon.style.marginTop = '0px';
  if(running) {
    resultDiv.innerHTML = '';
    resultDiv.append(icon);
  } else {
    resultDiv.innerHTML = '';
  }
}

function updateIcon (passed: boolean, iconDiv: Element) {
  const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
  icon.style.fontWeight = 'bold';
  icon.style.paddingRight = '5px';
  icon.style.marginTop = '0px';
  icon.style.color = passed ? 'lightgreen' : 'red';
  iconDiv.innerHTML = '';
  iconDiv.append(icon);
}

async function runTest(t: IPackageTest): Promise<boolean> {
  if (debugMode)
    debugger;
  const start = Date.now();
  testInProgress(t.resultDiv, true);
  const res = await grok.functions.call(
    `${t.packageName}:test`, {
    'category': t.test.category,
    'test': t.test.name,
  });
  const time = Date.now() - start;
  if (!testsResultsDf) {
    testsResultsDf = res;
    addPackageAndTimeInfo(testsResultsDf, time, t.packageName);
  } else {
    addPackageAndTimeInfo(res, time, t.packageName);
    removeTestRow(t.packageName, t.test.category, t.test.name);
    testsResultsDf = testsResultsDf.append(res);
  }
  const testRes = res.get('success', 0);
  updateTestResultsIcon(t.resultDiv, testRes);
  return testRes;
}

async function runAllTests(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
  testManagerView.path = '/' + testManagerView.name.replace(' ', '');
  switch (nodeType) {
    case NODE_TYPE.PACKAGE: {
      const progressBar = DG.TaskBarProgressIndicator.create(tests.package.name);
      testManagerView.path = `/${testManagerView.name.replace(' ', '')}/${tests.package.name}`;
      let testsSucceded = true;
      testInProgress(tests.resultDiv, true);
      await collectPackageTests(node as DG.TreeViewGroup, tests);
      const packageTests = packagesTests.filter(pt => pt.name === tests.package.name)[0];
      const cats = packageTests.categories;
      let completedTestsNum = 0;
      for (let cat of Object.values(cats)){
        testInProgress(tests.resultDiv, true);
        const res = await runTestsRecursive(cat, progressBar, packageTests.totalTests, completedTestsNum, tests.package.name);
        completedTestsNum = res.completedTests;
        if (!res.catRes)
          testsSucceded = false;
      }
      updateTestResultsIcon(tests.resultDiv, testsSucceded);
      progressBar.close();
      break;
    }
    case NODE_TYPE.CATEGORY: {
      const progressBar = DG.TaskBarProgressIndicator.create(`${tests.packageName}/${tests.fullName}`);
      testManagerView.path = `/${testManagerView.name.replace(' ', '')}/${tests.packageName}/${tests.fullName}`;
      await runTestsRecursive(tests, progressBar, tests.totalTests, 0, `${tests.packageName}/${tests.fullName}`);
      progressBar.close();
      break;
    }
    case NODE_TYPE.TEST: {
      testManagerView.path = `/${testManagerView.name.replace(' ', '')}/${tests.packageName}/${tests.test.category}/${tests.test.name}`;
      await runTest(tests);
      break;
    }
  }
  grok.shell.closeAll();
  grok.shell.v = testManagerView;
  setTimeout(() => {
    grok.shell.o = getTestsInfoPanel(node, tests, nodeType);
  }, 100);

}

async function runTestsRecursive(
  category: ICategory, 
  progressBar: DG.TaskBarProgressIndicator, 
  totalNumTests: number,
  completedTestsNum: number,
  progressInfo: string): Promise<any> {
  let testsSucceded = true;
  testInProgress(category.resultDiv, true);
  const subcats = Object.keys(category.subcategories);
  if (subcats.length > 0) {
    for (let subcat of subcats) {
      testInProgress(category.subcategories[subcat].resultDiv, true);
      const res = await runTestsRecursive(category.subcategories[subcat], progressBar, totalNumTests, completedTestsNum, progressInfo);
      if(!res.catRes)
        testsSucceded = false;
      completedTestsNum = res.completedTests;
    }
  }
  for (let t of category.tests) {
    const res = await runTest(t);
    if (!res)
    testsSucceded = false;
    completedTestsNum += 1;
    const percent = Math.floor(completedTestsNum/totalNumTests*100);
    progressBar.update(percent, `${progressInfo}: ${percent}% completed`);
  }
  updateTestResultsIcon(category.resultDiv, testsSucceded);
  return {completedTests: completedTestsNum, catRes: testsSucceded};
}

function addPackageAndTimeInfo (df: DG.DataFrame, time: number, pack: string) {
  df.columns.addNewInt('time, ms').init(() => time);
  df.columns.addNewString('package').init(() => pack);
}

function removeTestRow (pack: string, cat: string, test: string) {
  for (let i = 0; i < testsResultsDf.rowCount; i++) {
    if (testsResultsDf.get('package', i) === pack && testsResultsDf.get('category', i) === cat && testsResultsDf.get('name', i) === test) {
      testsResultsDf.rows.removeAt(i);
      return;
    }
  }
}

function getTestsInfoPanel(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
  const acc = ui.accordion();
  const accIcon = ui.element('i');
  accIcon.className = 'grok-icon svg-icon svg-view-layout';
  acc.addTitle(ui.span([accIcon, ui.label(`Tests details`)]));
  const grid = getTestsInfoGrid(resultsGridFilterCondition(tests, nodeType), nodeType);
  acc.addPane('Details', () => ui.div(testDetails(node, tests, nodeType)), true);
  acc.addPane('Results', () => ui.div(grid), true);
  return acc.root;
};

function resultsGridFilterCondition(tests: any, nodeType: NODE_TYPE) {
  return nodeType === NODE_TYPE.PACKAGE ? `Package = ${tests.package.name}` :
    nodeType === NODE_TYPE.CATEGORY ? `Package = ${tests.packageName} and category IN (${tests.subCatsNames.join(',')})` :
      `Package = ${tests.packageName} and category = ${tests.test.category} and name =  ${tests.test.name}`
}

function testDetails(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
  const detailsMap = nodeType === NODE_TYPE.PACKAGE ? { package: tests.package.name } :
    nodeType === NODE_TYPE.CATEGORY ? { package: tests.packageName, category: tests.fullName } :
      { package: tests.packageName, category: tests.test.category, test: tests.test.name };
  const detailsTable = ui.tableFromMap(detailsMap);
  const runButton = ui.bigButton('RUN', async () => {
    runAllTests(node, tests, nodeType);
  });
  return ui.divV([
    detailsTable,
    runButton
  ])
}


function getTestsInfoGrid(condition: string, nodeType: NODE_TYPE) {
  let info = ui.divText('No tests have been run');
  if (testsResultsDf) {
    const testInfo = testsResultsDf
      .groupBy(testsResultsDf.columns.names())
      .where(condition)
      .aggregate();
    if (testInfo.rowCount === 1 && testInfo.col('result').isNone(0))
      return info;
    if (nodeType === NODE_TYPE.TEST || testInfo.rowCount === 1 && !testInfo.col('result').isNone(0)) {
      const time = testInfo.get('time, ms', 0);
      const result = testInfo.get('result', 0);
      const resColor = testInfo.get('success', 0) ? 'lightgreen' : 'red';
      info = ui.divV([
        ui.divText(result, {style: {color: resColor}}),
        ui.divText(`Time, ms: ${time}`),
      ]);
      if (nodeType !== NODE_TYPE.TEST)
        info.append(ui.divText(`Test: ${testInfo.get('name', 0)}`));
        if (nodeType === NODE_TYPE.PACKAGE)
          info.append(ui.divText(`Category: ${testInfo.get('category', 0)}`));     
    } else
      info = ui.divV([
        ui.button('Add to workspace', () => {
          grok.shell.addTableView(testInfo);
        }),
        testInfo.plot.grid().root
      ])
  }
  return info;
};