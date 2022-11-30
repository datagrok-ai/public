import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {delay, Test, TestContext} from '@datagrok-libraries/utils/src/test';
import {c, testFunc} from './package';
import {Menu} from 'datagrok-api/dg';
import '../css/styles.css';

interface ITestManagerUI {
  testsTree: DG.TreeViewNode;
  runButton: HTMLButtonElement;
  runAllButton: HTMLButtonElement;
  debugButton: DG.InputBase<boolean>;
  benchmarkButton: DG.InputBase<boolean>;
  ribbonPanelDiv: HTMLDivElement;
}

interface IPackageTests {
  name: string;
  friendlyName: string;
  categories: {[index: string]: ICategory} | null;
  totalTests: number | null;
  resultDiv: HTMLElement;
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
  func: DG.Func | null
}

interface ITestFromUrl {
  packName: string,
  catName: string,
  testName: string
}

export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  grok.shell.addView(view);
  return view;
}

enum NODE_TYPE {
  PACKAGE = 'Package',
  CATEGORY = 'Category',
  TEST = 'TEST'
}

export class TestManager extends DG.ViewBase {
  packagesTests: IPackageTests[] = [];
  testsResultsDf: DG.DataFrame;
  testFunctions: any[];
  testManagerView: DG.View;
  selectedNode: DG.TreeViewGroup | DG.TreeViewNode;
  nodeDict: { [id: string]: any } = {};
  debugMode = false;
  benchmarkMode = false;
  tree: DG.TreeViewGroup;
  autoTestsCatName = 'Auto Tests';
  ribbonPanelDiv = undefined;
  dockLeft;

  constructor(name, dockLeft?: boolean) {
    super({});
    this.name = name;
    this.dockLeft = dockLeft;
  }

  async init(): Promise<void> {
    let pathSegments = window.location.pathname.split('/');
    pathSegments = pathSegments.map((it) => it ? it.replace(/%20/g, ' ') : undefined);
    this.testFunctions = await this.collectPackages();
    this.testManagerView = DG.View.create();
    const testFromUrl = pathSegments.length > 4 ?
      {packName: pathSegments[4], catName: pathSegments[5], testName: pathSegments[6]} : null;
    const testUIElements: ITestManagerUI = await this.createTestManagerUI(testFromUrl);
    this.testManagerView.name = this.name;
    addView(this.testManagerView);
    this.testManagerView.temp['ignoreCloseAll'] = true;
    this.ribbonPanelDiv = testUIElements.ribbonPanelDiv;
    this.testManagerView.append(testUIElements.ribbonPanelDiv);
    this.testManagerView.append(testUIElements.testsTree.root);
    if (this.dockLeft)
      grok.shell.dockManager.dock(this.testManagerView.root, 'left', null, this.name);
    this.runTestsForSelectedNode();
  }


  async collectPackages(packageName?: string): Promise<any[]> {
    let testFunctions = DG.Func.find({name: 'Test', meta: {file: 'package-test.js'}});
    testFunctions = testFunctions.sort((a, b) => a.package.friendlyName.localeCompare(b.package.friendlyName));
    if (packageName) testFunctions = testFunctions.filter((f: DG.Func) => f.package.name === packageName);
    return testFunctions;
  }

  async collectPackageTests(packageNode: DG.TreeViewGroup, f: any, testFromUrl?: ITestFromUrl) {
    if (this.testFunctions.filter((it) => it.package.name === f.package.name).length !== 0 &&
    this.packagesTests.filter((pt) => pt.name === f.package.name)[0].categories === null) {
      await f.package.load({file: f.options.file});
      const testModule = f.package.getModule(f.options.file);
      if (!testModule)
        console.error(`Error getting tests from '${f.package.name}/${f.options.file}' module.`);

      const allPackageTests = testModule ? testModule.tests : undefined;
      this.collectPackageAutoTests(f.package.name, allPackageTests);
      const packageTestsFinal: { [cat: string]: ICategory } = {};
      if (allPackageTests) {
        Object.keys(allPackageTests).forEach((cat) => {
          const tests: IPackageTest[] = allPackageTests[cat].tests.map((t) => {
            if (cat === this.autoTestsCatName) {
              const func = DG.Func.find({package: f.package.name, name: t.name})[0];
              return {test: t, packageName: f.package.name, func: func};
            } else
              return {test: t, packageName: f.package.name};
          });
          const subcats = cat.split(':');
          let subcatsFromUrl = [];
          if (testFromUrl && testFromUrl.catName)
            subcatsFromUrl = testFromUrl.catName.split(':');

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
            } else
              previousCat[subcats[i]].subCatsNames.push(cat);
            ;
            if (i === subcats.length - 1) {
              previousCat[subcats[i]].tests = previousCat[subcats[i]].tests ?
                previousCat[subcats[i]].tests.concat(tests) : tests;
            } else
              previousCat = previousCat[subcats[i]].subcategories;
          }
        });
        let testsNumInPack = 0;
        Object.keys(packageTestsFinal).forEach((cat) => {
          testsNumInPack += this.addCategoryRecursive(packageNode, packageTestsFinal[cat], testFromUrl);
        });
        const selectedPackage = this.packagesTests.filter((pt) => pt.name === f.package.name)[0];
        selectedPackage.categories = packageTestsFinal;
        selectedPackage.totalTests = testsNumInPack;
      };
    }
  }

  collectPackageAutoTests(packName: string, allPackageTests: any) {
    const packFunctions = DG.Func.find({package: packName});
    const autoTestsArray = [];
    for (const f of packFunctions) {
      const funcAutoTests = f.options['test'];
      if (funcAutoTests && Array.isArray(funcAutoTests) && funcAutoTests.length) {
        autoTestsArray.push(
          new Test(this.autoTestsCatName, f.name, null));
      }
    };
    if (autoTestsArray.length) {
      allPackageTests[this.autoTestsCatName] = {};
      allPackageTests[this.autoTestsCatName].tests = autoTestsArray;
    }
  }

  addCategoryRecursive(node: DG.TreeViewGroup, category: ICategory, testFromUrl?: ITestFromUrl) {
    const testPassed = ui.div();
    category.resultDiv = testPassed;
    const subnode = node.group(category.name, null, category.expand);
    subnode.root.children[0].append(testPassed);
    this.setRunTestsMenuAndLabelClick(subnode, category, NODE_TYPE.CATEGORY);
    if (testFromUrl && testFromUrl.catName === category.fullName)
      this.selectedNode = subnode;

    const subcats = Object.keys(category.subcategories);
    if (subcats.length > 0) {
      subcats.forEach((subcat) => {
        category.totalTests += this.addCategoryRecursive(subnode, category.subcategories[subcat], testFromUrl);
      });
    }
    category.tests.forEach((t) => {
      const testPassed = ui.div();
      const itemDiv = ui.divH([
        ui.divText(t.test.name),
        testPassed,
      ]);
      const item = subnode.item(itemDiv);
      t.resultDiv = testPassed;
      this.setRunTestsMenuAndLabelClick(item, t, NODE_TYPE.TEST);
      ui.tooltip.bind(item.root,
        () => this.getTestsInfoGrid(
          `Package = ${t.packageName} and category = ${t.test.category} and name =  ${t.test.name}`,
          NODE_TYPE.TEST, true));
      if (testFromUrl && testFromUrl.catName === category.fullName && testFromUrl.testName === t.test.name)
        this.selectedNode = item;
    });
    category.totalTests += category.tests.length;
    return category.totalTests;
  }


  async createTestManagerUI(testFromUrl: ITestFromUrl): Promise<ITestManagerUI> {
    this.tree = ui.tree();
    this.tree.root.classList.add('test-manager');
    this.tree.onSelectedNodeChanged.subscribe((res) => {
      this.selectedNode = res;
    });
    this.tree.root.style.width = '100%';
    this.tree.root.addEventListener('keyup', async (e) => {
      if (e.key === 'Enter')
        this.runTestsForSelectedNode();
    });
    for (const pack of this.testFunctions) {
      const testPassed = ui.div();
      this.packagesTests.push({
        name: pack.package.name,
        friendlyName: pack.package.friendlyName,
        categories: null,
        totalTests: null,
        resultDiv: testPassed,
      });
      const packNode = this.tree.group(pack.package.friendlyName,
        null, testFromUrl && pack.package.name === testFromUrl.packName);
      this.setRunTestsMenuAndLabelClick(packNode, pack, NODE_TYPE.PACKAGE);
      if (testFromUrl && testFromUrl.packName === pack.package.name) {
        this.selectedNode = packNode;
        await this.collectPackageTests(packNode, pack, testFromUrl);
      }
      packNode.root.children[0].append(testPassed);
      packNode.onNodeExpanding.subscribe(() => {
        this.collectPackageTests(packNode, pack);
      });
    }

    const {runAll, run, debug, benchmark} = this.createButtons();
    const {runAll: runAll1, run: run1, debug: debug1, benchmark: benchmark1} = this.createButtons();

    const ribbonPanelDiv = ui.divH([runAll1, run1, debug1.root, benchmark1.root],
      {style: {minHeight: '50px', maxHeight: '50px', alignItems: 'Center', borderBottom: '1px solid var(--grey-1)', paddingLeft: '5px'}});
    ribbonPanelDiv.classList.add('test');

    return {runAllButton: runAll, runButton: run, testsTree: this.tree,
      debugButton: debug, benchmarkButton: benchmark, ribbonPanelDiv: ribbonPanelDiv};
  }

  createButtons() {
    const runTestsButton = ui.button('Run', async () => {
      this.runTestsForSelectedNode();
    }, 'Run selected');

    const runAllButton = ui.bigButton('Run All', async () => {
      const nodes = this.tree.items;
      for (const node of nodes) {
        this.selectedNode = node;
        await this.runTestsForSelectedNode();
      }
    });

    //runAllButton.classList.add('btn-outline');
    runTestsButton.classList.add('ui-btn-outline');

    const debugButton = ui.boolInput('Debug', false, () => {this.debugMode = !this.debugMode;});
    debugButton.captionLabel.style.order = '1';
    debugButton.captionLabel.style.marginLeft = '5px';
    debugButton.root.style.marginLeft = '10px';

    const benchmarkButton = ui.boolInput('Benchmark', false, () => {this.benchmarkMode = !this.benchmarkMode;});
    benchmarkButton.captionLabel.style.order = '1';
    benchmarkButton.captionLabel.style.marginLeft = '5px';
    benchmarkButton.root.style.marginLeft = '5px';

    return {runAll: runAllButton, run: runTestsButton, debug: debugButton, benchmark: benchmarkButton};
  }

  async runTestsForSelectedNode() {
    if (this.selectedNode) {
      const id = this.selectedNode.root.id;
      await this.runAllTests(this.selectedNode, this.nodeDict[id].tests, this.nodeDict[id].nodeType);
    }
  }

  setRunTestsMenuAndLabelClick(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
    node.captionLabel.addEventListener('contextmenu', (e) => {
      Menu.popup()
        .item('Run', async () => {
          this.runAllTests(node, tests, nodeType);
        })
        .show();
      e.preventDefault();
      e.stopPropagation();
    });
    node.captionLabel.addEventListener('click', () => {
      grok.shell.o = this.getTestsInfoPanel(node, tests, nodeType);
    });
    if (nodeType === NODE_TYPE.TEST) {
      node.captionLabel.addEventListener('dblclick', async () => {
        this.runAllTests(node, tests, nodeType);
      });
    }
    const nodeId = Object.keys(this.nodeDict).length.toString();
    node.root.id = nodeId;
    this.nodeDict[Object.keys(this.nodeDict).length] = {tests: tests, nodeType: nodeType};
  }


  updateTestResultsIcon(resultDiv: HTMLElement, success?: boolean) {
    success === undefined ? resultDiv.innerHTML = '' : this.updateIcon(success, resultDiv);
  }

  testInProgress(resultDiv: HTMLElement, running: boolean) {
    const icon = ui.iconFA('spinner-third');
    icon.classList.add('fa-spin');
    icon.style.marginLeft = '2px';
    icon.style.marginTop = '0px';
    if (running) {
      resultDiv.innerHTML = '';
      resultDiv.append(icon);
    } else
      resultDiv.innerHTML = '';
  }

  updateIcon(passed: boolean, iconDiv: Element) {
    const icon = passed ? ui.iconFA('check') : ui.iconFA('times');
    icon.style.fontWeight = '500';
    icon.style.paddingLeft = '2px';
    icon.style.marginTop = '2px';
    icon.style.color = passed ? 'var(--green-2)' : 'var(--red-3)';
    iconDiv.innerHTML = '';
    iconDiv.append(icon);
  }

  async runTest(t: IPackageTest): Promise<boolean> {
    if (this.debugMode)
      debugger;
    this.testInProgress(t.resultDiv, true);
    const start = Date.now();
    let res;
    let testSucceeded;
    if (t.test.category === this.autoTestsCatName) {
      res = await testFunc(t.func);
      testSucceeded = Object.values(res as {[key: string]: {output: any, time: number}}).every((t) => t.output);
      res = this.createResultsDfForAutotests(res, t.func.name);
    } else {
      res = await grok.functions.call(
        `${t.packageName}:test`, {
          'category': t.test.category,
          'test': t.test.name,
          'testContext': new TestContext(false),
        });
      res.columns.addNewString('funcTest').init((i) => '');
      testSucceeded = res.get('success', 0);
    }
    const time = Date.now() - start;
    if (!this.testsResultsDf) {
      this.testsResultsDf = res;
      if (t.test.category !== this.autoTestsCatName)
        this.addTimeInfo(this.testsResultsDf, time);
      this.addPackageInfo(this.testsResultsDf, t.packageName);
    } else {
      if (t.test.category !== this.autoTestsCatName)
        this.addTimeInfo(res, time);
      if (res.col('package') == null)
        this.addPackageInfo(res, t.packageName);
      this.removeTestRow(t.packageName, t.test.category, t.test.name);
      this.testsResultsDf = this.testsResultsDf.append(res);
    }
    this.updateTestResultsIcon(t.resultDiv, testSucceeded);
    return testSucceeded;
  }

  createResultsDfForAutotests(res: {[key: string]: {output: any, time: number}}, funcName: string) {
    const resDf = DG.DataFrame.create(Object.keys(res).length);
    resDf.columns.addNewBool('success').init((i) => Object.values(res)[i].output === true);
    resDf.columns.addNewString('result').init((i) => Object.values(res)[i].output === true ? 'OK' : 'Failed');
    resDf.columns.addNewInt('ms').init((i) => Object.values(res)[i].time);
    resDf.columns.addNewString('category').init((i) => this.autoTestsCatName);
    resDf.columns.addNewString('name').init((i) => funcName);
    resDf.columns.addNewString('funcTest').init((i) => Object.keys(res)[i]);
    return resDf;
  }

  async runAllTests(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
    this.testManagerView.path = '/' + this.testManagerView.name.replace(' ', '');
    switch (nodeType) {
    case NODE_TYPE.PACKAGE: {
      const progressBar = DG.TaskBarProgressIndicator.create(tests.package.name);
      let testsSucceded = true;
      const packageTests = this.packagesTests.filter((pt) => pt.name === tests.package.name)[0];
      this.testInProgress(packageTests.resultDiv, true);
      await this.collectPackageTests(node as DG.TreeViewGroup, tests);
      const cats = packageTests.categories;
      let completedTestsNum = 0;
      for (const cat of Object.values(cats)) {
        this.testInProgress(packageTests.resultDiv, true);
        const res = await this.runTestsRecursive(cat, progressBar,
          packageTests.totalTests, completedTestsNum, tests.package.name);
        completedTestsNum = res.completedTests;
        if (!res.catRes)
          testsSucceded = false;
      }
      this.updateTestResultsIcon(packageTests.resultDiv, testsSucceded);
      this.testManagerView.path = `/${this.testManagerView.name.replace(' ', '')}/${tests.package.name}`;
      progressBar.close();
      break;
    }
    case NODE_TYPE.CATEGORY: {
      const progressBar = DG.TaskBarProgressIndicator.create(`${tests.packageName}/${tests.fullName}`);
      await this.runTestsRecursive(tests, progressBar, tests.totalTests, 0, `${tests.packageName}/${tests.fullName}`);
      this.testManagerView.path =
        `/${this.testManagerView.name.replace(' ', '')}/${tests.packageName}/${tests.fullName}`;
      progressBar.close();
      break;
    }
    case NODE_TYPE.TEST: {
      await this.runTest(tests);
      this.testManagerView.path =
        `/${this.testManagerView.name.replace(' ', '')}/${tests.packageName}/${tests.test.category}/${tests.test.name}`;
      break;
    }
    }
    await delay(1000);
    if (grok.shell.lastError.length > 0)
      grok.shell.error(`Unhandled exception: ${grok.shell.lastError}`);

    grok.shell.closeAll();
    setTimeout(() => {
      grok.shell.o = this.getTestsInfoPanel(node, tests, nodeType);
    }, 30);
  }

  async runTestsRecursive(
    category: ICategory,
    progressBar: DG.TaskBarProgressIndicator,
    totalNumTests: number,
    completedTestsNum: number,
    progressInfo: string): Promise<any> {
    let testsSucceded = true;
    this.testInProgress(category.resultDiv, true);
    const subcats = Object.keys(category.subcategories);
    if (subcats.length > 0) {
      for (const subcat of subcats) {
        this.testInProgress(category.subcategories[subcat].resultDiv, true);
        const res = await this.runTestsRecursive(
          category.subcategories[subcat], progressBar, totalNumTests, completedTestsNum, progressInfo);
        if (!res.catRes)
          testsSucceded = false;
        completedTestsNum = res.completedTests;
      }
    }
    for (const t of category.tests) {
      const res = await this.runTest(t);
      if (!res)
        testsSucceded = false;
      completedTestsNum += 1;
      const percent = Math.floor(completedTestsNum/totalNumTests*100);
      progressBar.update(percent, `${progressInfo}: ${percent}% completed`);
    }
    this.updateTestResultsIcon(category.resultDiv, testsSucceded);
    return {completedTests: completedTestsNum, catRes: testsSucceded};
  }

  addPackageInfo(df: DG.DataFrame, pack: string) {
    df.columns.addNewString('package').init(() => pack);
  }

  addTimeInfo(df: DG.DataFrame, time: number) {
    df.columns.addNewInt('time, ms').init(() => time);
  }

  removeTestRow(pack: string, cat: string, test: string) {
    for (let i = 0; i < this.testsResultsDf.rowCount; i++) {
      if (this.testsResultsDf.get('package', i) === pack &&
        this.testsResultsDf.get('category', i) === cat && this.testsResultsDf.get('name', i) === test) {
        this.testsResultsDf.rows.removeAt(i);
        return;
      }
    }
  }

  getTestsInfoPanel(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
    const acc = ui.accordion();
    const accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(`Tests details`)]));
    const grid = this.getTestsInfoGrid(this.resultsGridFilterCondition(tests, nodeType), nodeType);
    acc.addPane('Details', () => ui.div(this.testDetails(node, tests, nodeType)), true);
    acc.addPane('Results', () => ui.div(grid), true);
    return acc.root;
  };

  resultsGridFilterCondition(tests: any, nodeType: NODE_TYPE) {
    return nodeType === NODE_TYPE.PACKAGE ? `Package = ${tests.package.name}` :
      nodeType === NODE_TYPE.CATEGORY ?
        `Package = ${tests.packageName} and category IN (${tests.subCatsNames.join(',')})` :
        `Package = ${tests.packageName} and category = ${tests.test.category} and name = ${tests.test.name}`;
  }

  testDetails(node: DG.TreeViewGroup | DG.TreeViewNode, tests: any, nodeType: NODE_TYPE) {
    const detailsMap = nodeType === NODE_TYPE.PACKAGE ? {package: tests.package.name} :
      nodeType === NODE_TYPE.CATEGORY ? {package: tests.packageName, category: tests.fullName} :
        {package: tests.packageName, category: tests.test.category, test: tests.test.name};
    const detailsTable = ui.tableFromMap(detailsMap);
    const runButton = ui.bigButton('RUN', async () => {
      this.runAllTests(node, tests, nodeType);
    });
    return ui.divV([
      detailsTable,
      runButton,
    ]);
  }


  getTestsInfoGrid(condition: string, nodeType: NODE_TYPE, isTooltip?: boolean) {
    let info = ui.divText('No tests have been run');
    if (this.testsResultsDf) {
      const testInfo = this.testsResultsDf
        .groupBy(this.testsResultsDf.columns.names())
        .where(condition)
        .aggregate();
      if (testInfo.rowCount === 1 && testInfo.col('name').isNone(0))
        return info;
      const cat = testInfo.get('category', 0);
      if (testInfo.rowCount === 1 && !testInfo.col('name').isNone(0)) {
        const time = testInfo.get('time, ms', 0);
        const result = testInfo.get('result', 0);
        const resColor = testInfo.get('success', 0) ? 'var(--green-2)' : 'var(--red-3)';
        info = ui.divV([
          ui.divText(result, {style: {color: resColor, userSelect: 'text'}}),
          ui.divText(`Time, ms: ${time}`),
        ]);
        if (cat === this.autoTestsCatName)
          info.append(`Function: ${testInfo.get('funcTest', 0)}`);
        if (nodeType !== NODE_TYPE.TEST)
          info.append(ui.divText(`Test: ${testInfo.get('name', 0)}`));
        if (nodeType === NODE_TYPE.PACKAGE)
          info.append(ui.divText(`Category: ${cat}`));
      } else {
        if (!isTooltip) {
          info = ui.divV([
            ui.button('Add to workspace', () => {
              grok.shell.addTableView(testInfo);
            }),
            testInfo.plot.grid().root,
          ]);
        } else
          return null;
      }
    }
    return info;
  };
}


