import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EntityType} from './constants';
import {Test} from '@datagrok-libraries/utils/src/test';
import {c, _package} from './package';
import {delay} from '@datagrok-libraries/utils/src/test';

interface ITestManagerUI {
  runButton: HTMLButtonElement;
  testsTree: DG.TreeViewGroup;
}

interface IPackageTests {
  friendlyName: string;
  active: boolean;
  categories: {[index: string]: ICategory}
}

interface ICategory {
  active: boolean,
  tests: IPackageTest[]
}

interface IPackageTest {
  test: Test,
  active: boolean,
  packageName: string
}

interface ITestFromUrl {
  packName: string,
  catName: string,
  testName: string
}

let packagesTests: IPackageTests[];
let testsResultsDf: DG.DataFrame;
let testPath = '';

export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  view.path = '/' + view.name.replace(' ', '');
  grok.shell.addView(view);
  return view;
}

export async function testManagerView(): Promise<void> {
  let pathSegments = window.location.pathname.split('/');
  pathSegments = pathSegments.map((it) => it ? it.toLowerCase() : undefined);
  packagesTests = await collectTests({
    packName: pathSegments[4],
    catName: pathSegments[5],
    testName: pathSegments[6]});
  const v = DG.View.create();
  const testUIElements: ITestManagerUI = createTestManagerUI(true, v);
  v.name = 'Test Manager';
  addView(v);
  v.temp['ignoreCloseAll'] = true;
  v.setRibbonPanels(
    [[testUIElements.runButton]],
  );
  v.append(testUIElements.testsTree.root);
  runActiveTests(testUIElements.testsTree, v);
}


export async function _renderTestManagerPanel(ent: EntityType): Promise<DG.Widget> {
  if (ent == null)
    return DG.Widget.fromRoot(ui.divText('Entity does not exist.', {style: {color: 'var(--failure)'}}));

  if (ent.constructor.name === 'Package' || ent instanceof eval(`DG.Package`)) {
    packagesTests = await collectTests(undefined, ent.name);
    const testUIElements: ITestManagerUI = createTestManagerUI();
    const panelDiv = ui.divV([
      testUIElements.runButton,
      testUIElements.testsTree,
    ]);
    return DG.Widget.fromRoot(panelDiv);
  };
  return DG.Widget.fromRoot(ui.divText(`Entity is not of 'Package' type`, {style: {color: 'var(--failure)'}}));
}


async function collectTests(testFromUrl: ITestFromUrl, packageName?: string): Promise<IPackageTests[]> {
  let testFunctions = DG.Func.find({name: 'Test', meta: {file: 'package-test.js'}});
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  const packagesTestsList: IPackageTests[] = [];
  for (const f of testFunctions) {
    //await delay(4000); // will be removed with new version of js-api
    await f.package.load({file: f.options.file});
    const testModule = f.package.getModule(f.options.file);
    if (!testModule)
      console.error(`Error getting tests from '${f.package.name}/${f.options.file}' module.`);

    const allPackageTests = testModule ? testModule.tests : {};
    const testsWithPackNameAndActFlag: { [cat: string]: ICategory } = {};
    if (allPackageTests) {
      const normalizedPackName = f.package.name.toLowerCase();
      const packActive = testFromUrl ? testFromUrl.packName === normalizedPackName && !testFromUrl.catName : false;
      Object.keys(allPackageTests).forEach((cat) => {
        const catTests: IPackageTest[] = [];
        const normalizedCatName = cat.replace(/ /g, '').toLowerCase();
        const catActive = testFromUrl ?
          packActive || testFromUrl.packName === normalizedPackName && testFromUrl.catName === normalizedCatName && !testFromUrl.testName :
          false;
        allPackageTests[cat].tests
          .forEach((t) => catTests.push({
            test: t,
            active: testFromUrl ? catActive || testFromUrl.packName === normalizedPackName && testFromUrl.catName === normalizedCatName && testFromUrl.testName === t.name.replace(/ /g, '').toLowerCase() : false,
            packageName: f.package.name}));

        testsWithPackNameAndActFlag[cat] = {
          active: catActive,
          tests: catTests};
      });
      packagesTestsList.push({
        friendlyName: f.package.friendlyName,
        active: packActive,
        categories: testsWithPackNameAndActFlag});
    }
  }
  return packagesTestsList;
}


function createTestManagerUI(labelClick?: boolean, view?: DG.View): ITestManagerUI {
  const addCheckboxAndLabelClickListener = (item: DG.TreeViewNode, checked: boolean, isGroup: boolean, onChangeFunction: () => void, onItemClickFunction: () => void) => {
    item.enableCheckBox(checked);
    item.checkBox?.addEventListener('change', onChangeFunction);
    if (labelClick) {
      const label = isGroup ? item.root.children[0].children[2] : item.root.children[1];
      label.addEventListener('click', onItemClickFunction);
    }
  };

  const getTestsInfoAcc = (condition: string) => {
    const acc = ui.accordion();
    const accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(`Tests details`)]));
    const grid = getTestsInfoGrid(condition);
    acc.addPane('Results', () => ui.div(grid), true);
    return acc.root;
  };

  const getTestsInfoGrid = (condition: string) => {
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

  const tree = ui.tree();
  tree.root.style.width = '100%';
  packagesTests.forEach((pack) => {
    const packageGroup = tree.group(pack.friendlyName);
    addCheckboxAndLabelClickListener(packageGroup, pack.active, true, () => {
      pack.active = packageGroup.checked;
      Object.keys(pack.categories).forEach((cat) => {
        pack.categories[cat].active = packageGroup.checked;
        pack.categories[cat].tests.forEach((t) => t.active = packageGroup.checked);
      });
    },
    () => {
      grok.shell.o = getTestsInfoAcc(`Package = ${Object.values(pack.categories)[0].tests[0].packageName}`);
    });
    Object.keys(pack.categories).forEach((cat) => {
      const catGroup = packageGroup.group(cat);
      addCheckboxAndLabelClickListener(catGroup, pack.categories[cat].active, true, () => {
        pack.categories[cat].active = catGroup.checked;
        pack.categories[cat].tests.forEach((t) => {
          t.active = catGroup.checked;
        });
      },
      () => {
        grok.shell.o = getTestsInfoAcc(`Package = ${pack.categories[cat].tests[0].packageName} and category = ${cat}`);
      });
      pack.categories[cat].tests.forEach((t) => {
        const testPassed = ui.div();
        const itemDiv = ui.splitH([
          testPassed,
          ui.divText(t.test.name),
        ], {style: {display: 'block'}});
        const item = catGroup.item(itemDiv);
        item.root.id = `${t.packageName}|${cat}|${t.test.name}`;
        addCheckboxAndLabelClickListener(item, t.active, false, () => {
          t.active = item.checked;
        },
        () => {
          grok.shell.o = getTestsInfoAcc(`Package = ${t.packageName} and category = ${cat} and name =  ${t.test.name}`);
        });
        ui.tooltip.bind(item.root, () => getTestsInfoGrid(`Package = ${t.packageName} and category = ${cat} and name =  ${t.test.name}`));
      });
    });
  });

  const runTestsButton = ui.bigButton('Run', async () => {
    runActiveTests(tree, view);
  });

  return {runButton: runTestsButton, testsTree: tree};
}

function runActiveTests(tree: DG.TreeViewGroup, view?: DG.View) {
  const actTests = collectActiveTests();
  actTests.forEach((t) => updateTestResultsIcon(tree, t.packageName, t.test.category, t.test.name));
  if (actTests.length)
    runAllTests(actTests, tree, view);
}

function collectActiveTests() {
  let activeTests: IPackageTest[] = [];
  packagesTests.forEach((pack) => {
    Object.keys(pack.categories).forEach((cat) => {
      const active = pack.categories[cat].tests.filter((t) => t.active);
      if (active.length > 0) {

      }
      activeTests = activeTests.concat(pack.categories[cat].tests.filter((t) => t.active));
    });
  });
  const actPacks = [...new Set(activeTests.map((it) => it.packageName))];
  testPath = '';
  if (actPacks.length === 1) {
    testPath += `/${actPacks[0]}`;
    const actCats = [...new Set(activeTests.map((it) => it.test.category))];
    if (actCats.length === 1) {
      testPath += `/${actCats[0]}`;
      if (activeTests.length === 1)
        testPath += `/${activeTests[0].test.name}`;
    }
  }
  return activeTests;
}

function updateTestResultsIcon(tree: DG.TreeViewGroup, pack: string, cat: string, name: string, success?: boolean) {
  const items = tree.items;
  const item = items.filter((it) => it.root.id === `${pack}|${cat}|${name}`)[0];
  success === undefined ? item.root.children[1].children[0].children[0].innerHTML = '' : updateIcon(success, item.root.children[1].children[0].children[0]);
}

function updateIcon(passed: boolean, iconDiv: Element) {
  const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
  icon.style.fontWeight = 'bold';
  icon.style.paddingRight = '5px';
  icon.style.color = passed ? 'lightgreen' : 'red';
  iconDiv.append(icon);
}

async function runAllTests(activeTests: IPackageTest[], tree: DG.TreeViewGroup, view?: DG.View) {
  let completedTestsCount = 0;
  if (view) view.path = '/' + view.name.replace(' ', '') + testPath.replace(/ /g, '');
  for (let t of activeTests) {
    const start = Date.now();
    const res = await grok.functions.call(
      `${t.packageName}:test`, {
        'category': t.test.category,
        'test': t.test.name,
      });
      completedTestsCount +=1;
      if (!testsResultsDf) {
        testsResultsDf = res;
        addPackageAndTimeInfo(testsResultsDf, start, t.packageName);
      } else {
        addPackageAndTimeInfo(res, start, t.packageName);
        removeTestRow(t.packageName, t.test.category, t.test.name);
        testsResultsDf = testsResultsDf.append(res);
      }
      updateTestResultsIcon(tree, t.packageName, t.test.category, t.test.name, res.get('success', 0));
      if (completedTestsCount === activeTests.length) {
        grok.shell.closeAll();
        grok.shell.v = view;
      }
  }
};

function addPackageAndTimeInfo(df: DG.DataFrame, start: number, pack: string) {
  df.columns.addNewInt('time, ms').init(() => Date.now() - start);
  df.columns.addNewString('package').init(() => pack);
}

function removeTestRow(pack: string, cat: string, test: string) {
  for (let i = 0; i < testsResultsDf.rowCount; i++) {
    if (testsResultsDf.get('package', i) === pack && testsResultsDf.get('category', i) === cat && testsResultsDf.get('name', i) === test) {
      testsResultsDf.rows.removeAt(i);
      return;
    }
  }
}


