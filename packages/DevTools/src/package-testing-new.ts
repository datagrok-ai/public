import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EntityType } from './constants';
import { runTests, Test } from '@datagrok-libraries/utils/src/test';
import { c, _package } from './package';
import {delay} from "@datagrok-libraries/utils/src/test";
import { Menu } from 'datagrok-api/dg';

interface ITestManagerUI {
  testsTree: DG.TreeViewNode;
}

interface IPackageTests {
  name: string;
  friendlyName: string;
  active: boolean;
  categories: {[index: string]: ICategory}
}

interface ICategory {
  active: boolean,
  tests: IPackageTest[],
  subcategories: {[index: string]: ICategory},
  fullCatName: string
}

interface IPackageTest {
  test: Test,
  active: boolean,
  packageName: string,
  fullCatName: string,
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
  let v = DG.View.create();
  const testUIElements: ITestManagerUI = createTestManagerUI(true, v);
  v.name = 'Test Manager'
  addView(v);
  v.temp['ignoreCloseAll'] = true;
  v.append(testUIElements.testsTree.root); 
}

async function collectPackages(packageName?: string): Promise<any[]>  {
  let testFunctions = DG.Func.find({ name: 'Test' , meta: {file: 'package-test.js'}});
  if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
  return testFunctions;
}

async function collectPackageTests(packageNode: DG.TreeViewGroup, f: any) {
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
            previousCat[subcats[i]] = { active: true, tests: [], subcategories: {}, fullCatName: fullCatName };
          };
          if (i === subcats.length - 1) {
            previousCat[subcats[i]].tests = previousCat[subcats[i]].tests ? previousCat[subcats[i]].tests.concat(tests) : tests;
          } else {
            previousCat = previousCat[subcats[i]].subcategories;
          }
        }
      });
      Object.keys(packageTestsFinal).forEach(cat => {
        addCategoryRecursive(packageNode, packageTestsFinal[cat], cat);
      })
      packagesTests.push({
        name: f.package.name,
        friendlyName: f.package.friendlyName,
        active: true,
        categories: packageTestsFinal
      })
    };
  }
}

function addCategoryRecursive(node: DG.TreeViewGroup, category: ICategory, catName: string) {
  const subnode = node.group(catName, null, false);
  setRunTestsMenu(subnode, category, null);
  const subcats = Object.keys(category.subcategories);
  if (subcats.length > 0) {
    subcats.forEach((subcat) => {
      addCategoryRecursive(subnode, category.subcategories[subcat], subcat);
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
    setRunTestsMenu(item, null, t);
  })
}


function createTestManagerUI(labelClick?: boolean, view?: DG.View): ITestManagerUI {

  const tree = ui.tree();
  tree.root.style.width = '100%';
  testFunctions.forEach(pack => {
    const packNode = tree.group(pack.package.friendlyName, null, false);
    packNode.onNodeExpanding.subscribe(() => {
      collectPackageTests(packNode, pack);
    });
    setRunTestsMenu(packNode, null, null, pack);
  });


  return { testsTree: tree };
}

function setRunTestsMenu(node: DG.TreeViewGroup | DG.TreeViewNode, category: ICategory, test: IPackageTest, pack?: any) {
  node.captionLabel.addEventListener('contextmenu', (e) => {
    Menu.popup()
      .item('Run', async () => {      
        let testsToRun: IPackageTest[] = [];
        if(pack) {
          await collectPackageTests(node as  DG.TreeViewGroup, pack);
          const cats = packagesTests.filter(pt => pt.name === pack.package.name)[0].categories;
          Object.keys(cats).forEach(cat => {
            testsToRun = testsToRun.concat(collectTestsToRunRecursive(cats[cat]));
          })
        } else if (category){
          testsToRun = testsToRun.concat(collectTestsToRunRecursive(category));
        } else {
          testsToRun = [test];
        }
        runAllTests(testsToRun);
      })
      .show();
    e.preventDefault();
    e.stopPropagation();
  });
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

function updateIcon (passed: boolean, iconDiv: Element) {
  const icon = passed ? ui.iconFA('check') : ui.iconFA('ban');
  icon.style.fontWeight = 'bold';
  icon.style.paddingRight = '5px';
  icon.style.marginTop = '0px';
  icon.style.color = passed ? 'lightgreen' : 'red';
  iconDiv.innerHTML = '';
  iconDiv.append(icon);
}

async function runAllTests (activeTests: IPackageTest[], view?: DG.View) {
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
      updateTestResultsIcon(t.resultDiv, res.get('success', 0));
      if (completedTestsCount === activeTests.length) {
        grok.shell.closeAll();
       // grok.shell.v = view;
      }
  }
};

function addPackageAndTimeInfo (df: DG.DataFrame, start: number, pack: string) {
  df.columns.addNewInt('time, ms').init(() => Date.now() - start);
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