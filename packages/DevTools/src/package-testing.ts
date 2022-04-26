import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function _testManager() {
  let testFunctions = DG.Func.find({ name: 'Test' });
  const packagesTestsList = {};
  for (let f of testFunctions) {
    //@ts-ignore
    await f.package.load({ file: f.options.file });
    //@ts-ignore
    const allPackageTests = f.package.getModule(f.options.file).tests;
    if (allPackageTests) {
      //@ts-ignore
      packagesTestsList[ f.package.friendlyName ] = {name: f.package.name, tests: allPackageTests};
    }
  }
  createTestManagerUI(packagesTestsList);
}


function createTestManagerUI(packagesTests: any) {

  let testsResultsDf: DG.DataFrame;

  let addCheckboxAndLabelClickListener = (item: DG.TreeViewNode, isGroup: boolean, onChangeFunction: () => void, onItemClickFunction: () => void) => {
    item.enableCheckBox(false);
    item.checkBox?.addEventListener('change', onChangeFunction);
    const label = isGroup? item.root.children[0].children[2] : item.root.children[1];
    label.addEventListener('click', onItemClickFunction);
  };

  let applyToAllTests = (functionToApply: (t: any, p?: string) => void) => {
    Object.keys(packagesTests).forEach(pack => {
      Object.keys(packagesTests[ pack ].tests).forEach(cat => {
        //@ts-ignore
        packagesTests[ pack ].tests[ cat ].tests.forEach(t => functionToApply(t, packagesTests[ pack ].name));
      })
    });
  }

  applyToAllTests((t, p) => {
    t.packageName = p;
    t.active = false
  });

  let collectActiveTests = () => {
    //@ts-ignore
    let activeTests = [];
    Object.keys(packagesTests).forEach(pack => {
      Object.keys(packagesTests[ pack ].tests).forEach(cat => {
        //@ts-ignore
        activeTests = activeTests.concat(packagesTests[ pack ].tests[ cat ].tests.filter(t => t.active));
      });
    });
    //@ts-ignore
    return activeTests;
  }

  let addPackageAndTimeInfo = (df: DG.DataFrame, start: number, pack: string) => {
    df.columns.addNewInt('time, ms').init(() => Date.now() - start);
    df.columns.addNewString('package').init(() => pack);
  }

  let runAllTests = async (activeTests: any) => {
    //@ts-ignore
    activeTests.forEach(t => {
      const start = Date.now();
      grok.functions.call(
        `${t.packageName}:test`, {
        "category": t.category,
        "test": t.name
      }).then((res) => {
        if (!testsResultsDf) {
          testsResultsDf = res;
          addPackageAndTimeInfo(testsResultsDf, start, t.packageName);
        } else {
          addPackageAndTimeInfo(res, start, t.packageName);
          removeTestRow(t.packageName, t.category, t.name);
          testsResultsDf = testsResultsDf.append(res);
        }
        updateTestResultsIcon(tree, t.packageName, t.category, t.name, res.get('success', 0));
      })
    })
  };

  let removeTestRow = (pack, cat, test) => {
    for (let i = 0; i < testsResultsDf.rowCount; i++) {
      if (testsResultsDf.get('package', i) === pack && testsResultsDf.get('category', i) === cat && testsResultsDf.get('name', i) === test) {
        testsResultsDf.rows.removeAt(i);
        return;
      }
    }
  }


  let updateTestResultsIcon = (tree: DG.TreeViewNode, pack: string, cat: string, name: string, success?: boolean) => {
    const items = tree.items;
    const item = items.filter(it => it.root.id === `${pack}|${cat}|${name}`)[0];
    success === undefined ? item.root.children[1].children[0].children[0].innerHTML = '' : updateIcon(success, item.root.children[1].children[0].children[0]);
  }

  let updateIcon = (passed: boolean, iconDiv: Element) => {
    const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
    icon.style.fontWeight = 'bold';
    icon.style.paddingRight = '5px';
    icon.style.color = passed ? 'lightgreen' : 'red';
    iconDiv.append(icon);
  }

  let getTestsInfoDf = (condition: string) => {
    let acc = ui.accordion();
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(`Tests details`)]));
    let grid;
    if(testsResultsDf) {
      grid = testsResultsDf.
      groupBy(testsResultsDf.columns.names())
      .where(condition)
      .aggregate().plot.grid().root;
    }
    acc.addPane('Results', () => ui.div(grid), true);
    return acc.root;
  }

  const v = grok.shell.newView('Test manager');
  const tree = ui.tree();
  Object.keys(packagesTests).forEach(pack => {
    const packageGroup = tree.group(pack);
    addCheckboxAndLabelClickListener(packageGroup, true, () => {
      //@ts-ignore
      Object.keys(packagesTests[ pack ].tests).forEach(cat => {
        //@ts-ignore
        packagesTests[ pack ].tests[ cat ].tests.forEach(t => t.active = packageGroup.checked);
      });
    }, 
    () => {
      //@ts-ignore
      grok.shell.o = getTestsInfoDf(`Package = ${Object.values(packagesTests[ pack ].tests)[0].tests[0].packageName}`);
    });
    Object.keys(packagesTests[ pack ].tests).forEach(cat => {
      const catGroup = packageGroup.group(cat);
      addCheckboxAndLabelClickListener(catGroup, true, () => {
        packagesTests[ pack ].tests[ cat ].tests.forEach(t => {
          t.active = catGroup.checked;
        });
      },
      () => {
        grok.shell.o = getTestsInfoDf(`Package = ${packagesTests[ pack ].tests[cat].tests[0].packageName} and category = ${cat}`);
      });
      packagesTests[ pack ].tests[ cat ].tests.forEach(t => {
        let testPassed = ui.div();
        let itemDiv = ui.splitH([
          testPassed,
          ui.divText(t.name)
        ], false, {style: {display: 'block'}});
        let item = catGroup.item(itemDiv);
        item.root.id = `${t.packageName}|${cat}|${t.name}`;
        addCheckboxAndLabelClickListener(item, false, () => {
          t.active = item.checked;
        },
        () => {
          grok.shell.o = getTestsInfoDf(`Package = ${t.packageName} and category = ${cat} and name =  ${t.name}`);
        });
      });
    });
  });

  const runTestsButton = ui.bigButton('Run', async () => {
    let actTests = collectActiveTests();
    actTests.forEach(t => updateTestResultsIcon(tree, t.packageName, t.category, t.name));
    if(actTests.length) {
      runAllTests(actTests);
    }
  });

  v.setRibbonPanels(
    //@ts-ignore
    [ [ runTestsButton ] ] ,
  );
  v.append(tree.root);
}