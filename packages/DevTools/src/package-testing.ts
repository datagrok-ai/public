import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EntityType } from './constants';
import { Test } from '@datagrok-libraries/utils/src/test';

interface ITestManagerUI {
  runButton: HTMLButtonElement;
  testsTree: DG.TreeViewNode;
}

interface IPackageTests {
  friendlyName: string;
  categories: {[index: string]: IPackageTest[]}
}

interface IPackageTest {
  test: Test,
  active: boolean,
  packageName: string
}


export async function testManagerView(): Promise<void> {
  const packagesTestsList = await collectTests();
  const testUIElements: ITestManagerUI = createTestManagerUI(packagesTestsList, true);
  const v = grok.shell.newView('Test manager');
  v.setRibbonPanels(
    [ [ testUIElements.runButton ] ],
  );
  v.append(testUIElements.testsTree.root);
}


export async function _renderTestManagerPanel(ent: EntityType): Promise<DG.Widget> {
  if (ent == null) {
    return DG.Widget.fromRoot(ui.divText('Entity does not exist.', { style: { color: 'var(--failure)' } }));
  }
  if (ent.constructor.name === 'Package' || ent instanceof eval(`DG.Package`)) {
    const packagesTestsList = await collectTests(ent.name);
    const testUIElements: ITestManagerUI = createTestManagerUI(packagesTestsList);
    const panelDiv = ui.divV([
      testUIElements.runButton,
      testUIElements.testsTree
    ]);
    return DG.Widget.fromRoot(panelDiv);
  };
  return DG.Widget.fromRoot(ui.divText(`Entity is not of 'Package' type`, { style: { color: 'var(--failure)' } }));
}


async function collectTests(packageName?: string): Promise<IPackageTests[]>  {
  let testFunctions = DG.Func.find({ name: 'Test' });
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  const packagesTestsList: IPackageTests[] = [];
  for (let f of testFunctions) {
    await f.package.load({ file: f.options.file });
    const allPackageTests = f.package.getModule(f.options.file).tests;
    let testsWithPackNameAndActFlag: { [cat: string]: IPackageTest[] } = {};
    if (allPackageTests) {
      Object.keys(allPackageTests).forEach((cat) => {
        let catTests: IPackageTest[] = [];
        allPackageTests[cat].tests.forEach((t) => catTests.push({test: t, active: false, packageName: f.package.name}));
        testsWithPackNameAndActFlag[cat] = catTests;
      });
      packagesTestsList.push({friendlyName: f.package.friendlyName, categories: testsWithPackNameAndActFlag});
    }
  }
  return packagesTestsList;
}


function createTestManagerUI(packagesTests: IPackageTests[], labelClick?: boolean): ITestManagerUI {

  let testsResultsDf: DG.DataFrame;

  let addCheckboxAndLabelClickListener = (item: DG.TreeViewNode, isGroup: boolean, onChangeFunction: () => void, onItemClickFunction: () => void) => {
    item.enableCheckBox(false);
    item.checkBox?.addEventListener('change', onChangeFunction);
    if (labelClick) {
      const label = isGroup ? item.root.children[ 0 ].children[ 2 ] : item.root.children[ 1 ];
      label.addEventListener('click', onItemClickFunction);
    }
  };

  let collectActiveTests = () => {
    let activeTests: IPackageTest[]= [];
    packagesTests.forEach(pack => {
      Object.keys(pack.categories).forEach(cat => {
        activeTests = activeTests.concat(pack.categories[ cat ].filter(t => t.active));
      });
    });
    return activeTests;
  }

  let addPackageAndTimeInfo = (df: DG.DataFrame, start: number, pack: string) => {
    df.columns.addNewInt('time, ms').init(() => Date.now() - start);
    df.columns.addNewString('package').init(() => pack);
  }

  let runAllTests = async (activeTests: IPackageTest[]) => {
    activeTests.forEach(t => {
      const start = Date.now();
      grok.functions.call(
        `${t.packageName}:test`, {
        "category": t.test.category,
        "test": t.test.name
      }).then((res) => {
        if (!testsResultsDf) {
          testsResultsDf = res;
          addPackageAndTimeInfo(testsResultsDf, start, t.packageName);
        } else {
          addPackageAndTimeInfo(res, start, t.packageName);
          removeTestRow(t.packageName, t.test.category, t.test.name);
          testsResultsDf = testsResultsDf.append(res);
        }
        updateTestResultsIcon(tree, t.packageName, t.test.category, t.test.name, res.get('success', 0));
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
  packagesTests.forEach(pack => {
    const packageGroup = tree.group(pack.friendlyName);
    addCheckboxAndLabelClickListener(packageGroup, true, () => {
      Object.keys(pack.categories).forEach(cat => {
        pack.categories[cat].forEach(t => t.active = packageGroup.checked);
      });
    },
      () => {
        grok.shell.o = getTestsInfoAcc(`Package = ${Object.values(pack.categories)[0][0].packageName}`);
      });
    Object.keys(pack.categories).forEach(cat => {
      const catGroup = packageGroup.group(cat);
      addCheckboxAndLabelClickListener(catGroup, true, () => {
        pack.categories[ cat ].forEach(t => {
          t.active = catGroup.checked;
        });
      },
        () => {
          grok.shell.o = getTestsInfoAcc(`Package = ${pack.categories[ cat ][ 0 ].packageName} and category = ${cat}`);
        });
      pack.categories[ cat ].forEach(t => {
        let testPassed = ui.div();
        let itemDiv = ui.splitH([
          testPassed,
          ui.divText(t.test.name)
        ], { style: { display: 'block' } });
        let item = catGroup.item(itemDiv);
        item.root.id = `${t.packageName}|${cat}|${t.test.name}`;
        addCheckboxAndLabelClickListener(item, false, () => {
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
