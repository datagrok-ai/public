import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EntityType } from './constants';
import { runTests, Test } from '@datagrok-libraries/utils/src/test';
import { c, _package } from './package';
import {delay} from "@datagrok-libraries/utils/src/test";
import { Menu } from 'datagrok-api/dg';

enum NODE_TYPE {
  PACKAGE = 'Package',
  CATEGORY = 'Category',
  TEST = 'TEST'
}

interface ITestManagerUI {
  testsTree: DG.TreeViewNode;
}

interface IPackageTests {
  name: string;
  friendlyName: string;
  categories: {[index: string]: ICategory}
}

interface ICategory {
  name: string,
  tests: IPackageTest[],
  subcategories: {[index: string]: ICategory},
  packageName: string,
  subCatsNames: string[],
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
let testPath = '';
let testFunctions: any[];
let testManagerView: DG.View;
let selectedNode: DG.TreeViewGroup|DG.TreeViewNode;
const nodeDict: {[id: string]: any} = {};

export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  view.path = '/' + view.name.replace(' ', '');
  grok.shell.addView(view);
  return view;
}

export async function testManagerViewNew(): Promise<void> {
  let pathSegments = window.location.pathname.split('/');
  pathSegments = pathSegments.map(it => it ? it.toLowerCase() : undefined);
  testFunctions = await collectPackages();
  testManagerView = DG.View.create();
  const testUIElements: ITestManagerUI = createTestManagerUI();
  testManagerView.name = 'Test Manager'
  addView(testManagerView);
  testManagerView.temp['ignoreCloseAll'] = true;
  testManagerView.append(testUIElements.testsTree.root); 
}

async function collectPackages(packageName?: string): Promise<any[]>  {
  let testFunctions = DG.Func.find({ name: 'Test' , meta: {file: 'package-test.js'}});
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  return testFunctions;
}

async function collectPackageTests(packageNode: DG.TreeViewGroup, f: any, expand: boolean) {
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
        let previousCat = packageTestsFinal;
        for (let i = 0; i < subcats.length; i++) {
          if (!packageTestsFinal[subcats[i]]) {
            previousCat[subcats[i]] = { tests: [], subcategories: {}, packageName: f.package.name, name: subcats[i], subCatsNames: [cat] };
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
      Object.keys(packageTestsFinal).forEach(cat => {
        addCategoryRecursive(packageNode, packageTestsFinal[cat], expand);
      })
      packagesTests.push({
        name: f.package.name,
        friendlyName: f.package.friendlyName,
        categories: packageTestsFinal
      })
    };
  }
}

function addCategoryRecursive(node: DG.TreeViewGroup, category: ICategory, expand: boolean) {
  const subnode = node.group(category.name, null, expand);
  setRunTestsMenuAndLabelClick(subnode, category, NODE_TYPE.CATEGORY);
  const subcats = Object.keys(category.subcategories);
  if (subcats.length > 0) {
    subcats.forEach((subcat) => {
      addCategoryRecursive(subnode, category.subcategories[subcat], expand);
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
  })
}


function createTestManagerUI(): ITestManagerUI {
  const tree = ui.tree();
  tree.onSelectedNodeChanged.subscribe((res) => {
    selectedNode = res;
  });
  tree.root.style.width = '100%';
  tree.root.addEventListener('keyup', async (e) => {
    if(e.key === 'Enter') {
      if (selectedNode) {
      const id = selectedNode.root.id;
      //@ts-ignore
      let testsToRun = await collectTestsToRun(selectedNode, nodeDict[id].tests, nodeDict[id].nodeType);
        runAllTests(testsToRun);
      }
    }
  });
  testFunctions.forEach(pack => {
    const packNode = tree.group(pack.package.friendlyName, null, false);
    packNode.onNodeExpanding.subscribe(() => {
      collectPackageTests(packNode, pack, false);
    });
    setRunTestsMenuAndLabelClick(packNode, pack, NODE_TYPE.PACKAGE);
  });


  return { testsTree: tree };
}

function setRunTestsMenuAndLabelClick(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
  node.captionLabel.addEventListener('contextmenu', (e) => {
    Menu.popup()
      .item('Run', async () => {
        let testsToRun = await collectTestsToRun(node, tests, nodeType);
        runAllTests(testsToRun);
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
      let testsToRun = await collectTestsToRun(node, tests, nodeType);
          runAllTests(testsToRun);
    });
  }
    const nodeId = Object.keys(nodeDict).length.toString();
    node.root.id = nodeId;
    nodeDict[Object.keys(nodeDict).length] = {tests: tests, nodeType: nodeType};
}

async function collectTestsToRun(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE): Promise<IPackageTest[]> {
  let testsToRun: IPackageTest[] = [];
  switch (nodeType) {
    case NODE_TYPE.PACKAGE: {
      await collectPackageTests(node as DG.TreeViewGroup, tests, true);
      const cats = packagesTests.filter(pt => pt.name === tests.package.name)[0].categories;
      Object.keys(cats).forEach(cat => {
        testsToRun = testsToRun.concat(collectTestsToRunRecursive(cats[cat]));
      });
      break;
    }
    case NODE_TYPE.CATEGORY: {
      testsToRun = testsToRun.concat(collectTestsToRunRecursive(tests));
      break;
    }
    case NODE_TYPE.TEST: {
      testsToRun = [tests];
      break;
    }
  }
  return testsToRun;
}


function collectTestsToRunRecursive(category: ICategory) {
  const subcats = Object.keys(category.subcategories);
  if (subcats.length > 0) {
    let tests = [];
    subcats.forEach((subcat) => {
      tests = tests.concat(collectTestsToRunRecursive(category.subcategories[subcat]));
    })
    return tests.concat(category.tests);
  } 
  return category.tests;
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

async function runAllTests(activeTests: IPackageTest[], view?: DG.View) {
  let completedTestsCount = 0;
  if (view) view.path = '/' + view.name.replace(' ', '') + testPath.replace(/ /g, '');
  for (let t of activeTests) {
    const start = Date.now();
    testInProgress(t.resultDiv, true);
    const res = await grok.functions.call(
      `${t.packageName}:test`, {
      'category': t.test.category,
      'test': t.test.name,
    });
    const time = Date.now() - start;
    testInProgress(t.resultDiv, false);
    completedTestsCount += 1;
    if (!testsResultsDf) {
      testsResultsDf = res;
      addPackageAndTimeInfo(testsResultsDf, time, t.packageName);
    } else {
      addPackageAndTimeInfo(res, time, t.packageName);
      removeTestRow(t.packageName, t.test.category, t.test.name);
      testsResultsDf = testsResultsDf.append(res);
    }
    updateTestResultsIcon(t.resultDiv, res.get('success', 0));
    if (completedTestsCount === activeTests.length) {
      grok.shell.closeAll();
      grok.shell.v = testManagerView;
    }
  }
};

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
  const grid = getTestsInfoGrid(resultsGridFilterCondition(tests, nodeType));
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
  const detailsMap = nodeType === NODE_TYPE.PACKAGE ? { package: tests.package.name, categories: 'all', tests: 'all' } :
    nodeType === NODE_TYPE.CATEGORY ? { package: tests.packageName, categories: tests.subCatsNames.join(','), tests: 'all' } :
      { package: tests.packageName, category: tests.test.category, test: tests.test.name };
  const detailsTable = ui.tableFromMap(detailsMap);
  const runButton = ui.bigButton('RUN', async () => {
    let testsToRun = await collectTestsToRun(node, tests, nodeType);
    runAllTests(testsToRun);
  });
  return ui.divV([
    detailsTable,
    runButton
  ])
}


function getTestsInfoGrid(condition: string) {
  let grid;
  if (testsResultsDf) {
    const testInfo = testsResultsDf
      .groupBy(testsResultsDf.columns.names())
      .where(condition)
      .aggregate();
    if (testInfo.rowCount === 1 && !testInfo.col('result').isNone(0)) {
      const gridMap = {};
      testInfo.columns.names().forEach((col) => {
        gridMap[col] = testInfo.get(col, 0);
      });
      grid = ui.tableFromMap(gridMap);
    } else
      grid = testInfo.plot.grid().root;
  }
  return grid;
};