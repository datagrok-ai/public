import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EntityType } from './constants';
import { Test } from '@datagrok-libraries/utils/src/test';
import { c, _package } from './package';
import {delay} from "@datagrok-libraries/utils/src/test";

interface ITestManagerUI {
  runButton: HTMLButtonElement;
  testsTree: DG.TreeViewNode;
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


export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  view.path = '/' + view.name.replace(' ', '');
  grok.shell.addView(view);
  return view;
}

export async function testManagerView(): Promise<void> {
  let pathSegments = window.location.pathname.split('/');
  pathSegments = pathSegments.map(it => it ? it.toLowerCase() : undefined);
  const packagesTestsList = await collectTests({
    packName: pathSegments[4], 
    catName: pathSegments[5], 
    testName: pathSegments[6]});
  const testUIElements: ITestManagerUI = createTestManagerUI(packagesTestsList, true);
  let v = DG.View.create();
  v.name = 'Test Manager'
  addView(v);
  v.temp['ignoreCloseAll'] = true;
  v.setRibbonPanels(
    [ [ testUIElements.runButton ] ],
  );
  v.append(testUIElements.testsTree.root);
  testUIElements.runButton.click();  
}


export async function _renderTestManagerPanel(ent: EntityType): Promise<DG.Widget> {
  if (ent == null) {
    return DG.Widget.fromRoot(ui.divText('Entity does not exist.', { style: { color: 'var(--failure)' } }));
  }
  if (ent.constructor.name === 'Package' || ent instanceof eval(`DG.Package`)) {
    const packagesTestsList = await collectTests(undefined, ent.name);
    const testUIElements: ITestManagerUI = createTestManagerUI(packagesTestsList);
    const panelDiv = ui.divV([
      testUIElements.runButton,
      testUIElements.testsTree
    ]);
    return DG.Widget.fromRoot(panelDiv);
  };
  return DG.Widget.fromRoot(ui.divText(`Entity is not of 'Package' type`, { style: { color: 'var(--failure)' } }));
}


async function collectTests(testFromUrl: ITestFromUrl, packageName?: string): Promise<IPackageTests[]>  {
  let testFunctions = DG.Func.find({ name: 'Test' });
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  const packagesTestsList: IPackageTests[] = [];
  for (let f of testFunctions) {
    await f.package.load({ file: f.options.file });
   // await delay(2000);
    const allPackageTests = f.package.getModule(f.options.file).tests;
    let testsWithPackNameAndActFlag: { [cat: string]: ICategory } = {};
    if (allPackageTests) {
      const normalizedPackName = f.package.name.toLowerCase();
      const packActive = testFromUrl ? testFromUrl.packName === normalizedPackName && !testFromUrl.catName : false;
      Object.keys(allPackageTests).forEach((cat) => {
        let catTests: IPackageTest[] = [];
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


function createTestManagerUI(packagesTests: IPackageTests[], labelClick?: boolean): ITestManagerUI {

  let testsResultsDf: DG.DataFrame;

  let addCheckboxAndLabelClickListener = (item: DG.TreeViewNode, checked: boolean, isGroup: boolean, onChangeFunction: () => void, onItemClickFunction: () => void) => {
    item.enableCheckBox(checked);
    item.checkBox?.addEventListener('change', onChangeFunction);
    if (labelClick) {
      const label = isGroup ? item.root.children[ 0 ].children[ 2 ] : item.root.children[ 1 ];
      label.addEventListener('click', onItemClickFunction);
    }
  };

  let collectActiveTests = () => {
    let activeTests: IPackageTest[] = [];
    packagesTests.forEach(pack => {
      Object.keys(pack.categories).forEach(cat => {
        activeTests = activeTests.concat(pack.categories[ cat ].tests.filter(t => t.active));
      });
    });
    return activeTests;
  }

  let addPackageAndTimeInfo = (df: DG.DataFrame, start: number, pack: string) => {
    df.columns.addNewInt('time, ms').init(() => Date.now() - start);
    df.columns.addNewString('package').init(() => pack);
  }

  let runAllTests = async (activeTests: IPackageTest[]) => {
    let completedTestsCount = 0;
    activeTests.forEach(t => {
      const start = Date.now();
      grok.functions.call(
        `${t.packageName}:test`, {
        "category": t.test.category,
        "test": t.test.name
      }).then((res) => {
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
        if(completedTestsCount === activeTests.length) {
        //  grok.shell.closeAll();
        }
      })
    })
  };

  let removeTestRow = (pack: string, cat: string, test: string) => {
    for (let i = 0; i < testsResultsDf.rowCount; i++) {
      if (testsResultsDf.get('package', i) === pack && testsResultsDf.get('category', i) === cat && testsResultsDf.get('name', i) === test) {
        testsResultsDf.rows.removeAt(i);
        return;
      }
    }
  }

  let updateTestResultsIcon = (tree: DG.TreeViewNode, pack: string, cat: string, name: string, success?: boolean) => {
    const items = tree.items;
    const item = items.filter(it => it.root.id === `${pack}|${cat}|${name}`)[ 0 ];
    success === undefined ? item.root.children[ 1 ].children[ 0 ].children[ 0 ].innerHTML = '' : updateIcon(success, item.root.children[ 1 ].children[ 0 ].children[ 0 ]);
  }

  let updateIcon = (passed: boolean, iconDiv: Element) => {
    const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
    icon.style.fontWeight = 'bold';
    icon.style.paddingRight = '5px';
    icon.style.color = passed ? 'lightgreen' : 'red';
    iconDiv.append(icon);
  }

  let getTestsInfoAcc = (condition: string) => {
    let acc = ui.accordion();
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([ accIcon, ui.label(`Tests details`) ]));
    let grid = getTestsInfoGrid(condition);
    acc.addPane('Results', () => ui.div(grid), true);
    return acc.root;
  }

  let getTestsInfoGrid = (condition: string) => {
    let grid;
    if (testsResultsDf) {
      let testInfo = testsResultsDf
        .groupBy(testsResultsDf.columns.names())
        .where(condition)
        .aggregate();
      if (testInfo.rowCount === 1 && !testInfo.col('result').isNone(0)) {
        const gridMap = {};
        testInfo.columns.names().forEach(col => {
          gridMap[ col ] = testInfo.get(col, 0);
        });
        grid = ui.tableFromMap(gridMap);
      } else {
        grid = testInfo.plot.grid().root;
      }
    }
    return grid;
  }

  const tree = ui.tree();
  tree.root.style.width = '100%';
  packagesTests.forEach(pack => {
    const packageGroup = tree.group(pack.friendlyName);
    addCheckboxAndLabelClickListener(packageGroup, pack.active, true, () => {
      pack.active = packageGroup.checked;
      Object.keys(pack.categories).forEach(cat => {
        pack.categories[cat].active = packageGroup.checked;
        pack.categories[cat].tests.forEach(t => t.active = packageGroup.checked);
      });
    },
      () => {
        grok.shell.o = getTestsInfoAcc(`Package = ${Object.values(pack.categories)[0].tests[0].packageName}`);
      });
    Object.keys(pack.categories).forEach(cat => {
      const catGroup = packageGroup.group(cat);
      addCheckboxAndLabelClickListener(catGroup, pack.categories[cat].active, true, () => {
        pack.categories[cat].active = catGroup.checked;
        pack.categories[ cat ].tests.forEach(t => {
          t.active = catGroup.checked;
        });
      },
        () => {
          grok.shell.o = getTestsInfoAcc(`Package = ${pack.categories[ cat ].tests[ 0 ].packageName} and category = ${cat}`);
        });
      pack.categories[ cat ].tests.forEach(t => {
        let testPassed = ui.div();
        let itemDiv = ui.splitH([
          testPassed,
          ui.divText(t.test.name)
        ], { style: { display: 'block' } });
        let item = catGroup.item(itemDiv);
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
    let actTests = collectActiveTests();
    actTests.forEach(t => updateTestResultsIcon(tree, t.packageName, t.test.category, t.test.name));
    if (actTests.length) {
      runAllTests(actTests);
    }
  });

  return { runButton: runTestsButton, testsTree: tree };
}
