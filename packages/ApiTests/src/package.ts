/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/data_frame';
import './dataframe/calculated-columns';
import './dataframe/events';
import './shell/shell';
import './shell/windows';
import './viewer/viewer';
import './views/docking';
import './views/events';
import './views/layouts';
import './dapi/files';
import './dapi/fetch';
import './dapi/admin';
import './dapi/groups';
import './ui/inputs';
import './ui/forms';
import './dapi/dapi';
import './dapi/entities';
import './dapi/layouts';
import './dapi/projects';
import './dapi/tables';
import './dapi/user-data-storage';
import './dapi/users';
import './shell/ml';
import './ui/divs';
import './ui/buttons';
import './widgets/legend';
import './ui/icons';
import './ui/tables';
import './ui/rangeSlider';
import './ui/accordion';
import './ui/tabControl';
import './ui/list';
import './ui/image';
import './ui/viewers-adding'
import './grid/grid';

import {runTests} from "@datagrok-libraries/utils/src/test";
export let _package = new DG.Package();


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function test(category: string, test: string): Promise<DG.DataFrame> {
  let data = await runTests({category, test});
  return DG.DataFrame.fromObjects(data)!;
}

//name: testPackages
//output: dataframe result
//top-menu: Tools | Dev | Test Packages
export async function testPackages(): Promise<DG.DataFrame> {
  let funcs = DG.Func.find({name:'test'});
  let dfs:DG.DataFrame[] = [];
  for (const f of funcs) {
    if (f.package?.name != null) {
      grok.shell.closeAll();
      grok.shell.info(`Testing ${f.package.name}`);
      let df = await f.apply();
      if (df == null) {
        grok.shell.error(`Failed to fetch test results from ${f.package.name}`);
        continue;
      }
      let packageColumn = DG.Column.string('package', df.rowCount);
      packageColumn.init((n) => f.package.name);
      df.columns.insert(packageColumn, 0);
      dfs.push(df);
      grok.shell.closeAll();
    }
  }

  let result: DG.DataFrame | null = null;
  for (const df of dfs) {
    if (result == null)
      result = df;
    else result.append(df, true);
  }

  return result!;
}



//name: testManager
//top-menu: Tools | Dev | Test Manager
export async function testManager() {
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
  createTestMangerUI(packagesTestsList);
}


function createTestMangerUI(packagesTests: any) {

  let addCheckbox = (item: any, onChangeFunction: () => void) => {
    item.enableCheckBox(false);
    item.checkBox?.addEventListener('change', onChangeFunction);
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

  let runAllTests = async (activeTests: any) => {
    //@ts-ignore
    const promises = activeTests.map((t) => 
      grok.functions.call(
        `${t.packageName}:test`, {
        "category": t.category,
        "test": t.name
      })
    );
    Promise.all(promises).then((res) => {
      const testsResultsDf = res.reduce(
        (previousDf, currentDf) => previousDf.append(currentDf)
      );
      grok.shell.addTableView(testsResultsDf);
      updateTestResultsIcons(tree, testsResultsDf);
    });
  };

  let updateTestResultsIcons = (tree: DG.TreeViewNode, testResults: DG.DataFrame) => {
    const items = tree.items;
    let rowCount = testResults.rowCount;
    for (let i = 0; i < rowCount; i++){
        const item = items.filter(it => it.root.id === `${testResults.get('category', i)}|${testResults.get('name', i)}`)[0];
        updateIcon(testResults.get('success', i), item.root.children[1].children[0].children[1]);
    }
  }

  let updateIcon = (passed: boolean, iconDiv: Element) => {
    const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
    icon.style.fontWeight = 'bold';
    icon.style.paddingLeft = '5px';
    icon.style.color = passed ? 'lightgreen' : 'red';
    iconDiv.innerHTML = '';
    iconDiv.append(icon);
  }

  const v = grok.shell.newView('Test manager');
  const tree = ui.tree();
  Object.keys(packagesTests).forEach(pack => {
    const packageGroup = tree.group(pack);
    addCheckbox(packageGroup, () => {
      //@ts-ignore
      Object.keys(packagesTests[ pack ].tests).forEach(cat => {
        //@ts-ignore
        packagesTests[ pack ].tests[ cat ].tests.forEach(t => t.active = packageGroup.checked);
      });
    });
    Object.keys(packagesTests[ pack ].tests).forEach(cat => {
      const catGroup = packageGroup.group(cat);
      addCheckbox(catGroup, () => {
        //@ts-ignore
        packagesTests[ pack ].tests[ cat ].tests.forEach(t => {
          t.active = catGroup.checked;
        });
      });
      //@ts-ignore
      packagesTests[ pack ].tests[ cat ].tests.forEach(t => {
        let testPassed = ui.div();
        let itemDiv = ui.splitH([
          ui.divText(t.name),
          testPassed
        ], {style: {display: 'block'}});
        let item = catGroup.item(itemDiv);
        item.root.id = `${cat}|${t.name}`;
        addCheckbox(item, () => {
          t.active = item.checked;
        });
      });
    });
  });

  const runTestsButton = ui.bigButton('Run', async () => {
    let actTests = collectActiveTests();
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

