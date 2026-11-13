import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UsageWidget} from './widgets/usage-widget';
import {PackageUsageWidget} from './widgets/package-usage-widget';
import '../css/usage_analysis.css';
import '../css/test_track.css';
import '../css/click_events_widget.css';
import {ViewHandler} from './view-handler';
import {TestTrack} from './test-track/app';
import {ReportsWidget} from './widgets/reports-widget';
import {ReportingApp} from './reporting/reporting_app';
import {TestAnalysisManager} from './test-analysis/test-analysis-manager';
import {getDate} from './utils';
import {TestGridCellHandler} from './handlers/test-grid-cell-handler';
import {initTestStickyMeta} from './test-analysis/sticky-meta-initialization';
import {TestDashboardWidget} from './viewers/ua-test-dashboard-viewer';
import {ClickEventsWidget} from './widgets/click-events-widget';
import {fetchVexDashboardImages, vexImagesToDataFrame} from './tabs/vulnerabilities';
import {ReleaseHandler} from './release/release-handler';
import {isDevHost} from './release/data';

export const _package = new DG.Package();
export * from './package.g';


export class PackageFunctions {
  @grok.decorators.init()
  static _initUA(): void {
    DG.ObjectHandler.register(new TestGridCellHandler());
    initTestStickyMeta();
  }


  @grok.decorators.func({
    'meta': {
      'url': '/tests/list',
    },
  })
  static async TestsList(): Promise<DG.DataFrame| undefined> {
    const pacakageTests = await TestAnalysisManager.collectPackageTests();
    const packageTestsListMapped = pacakageTests.map((elem) => {
      return {'name': 'test-package ' + elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const manualTest = await TestAnalysisManager.collectManualTestNames();
    const manualTestsListMapped = manualTest.map((elem) => {
      return {'name': 'test-manual ' + elem};
    });
    const resultTestsList = manualTestsListMapped.concat(packageTestsListMapped);
    return DG.DataFrame.fromObjects(resultTestsList);
  }


  @grok.decorators.func({
    'meta': {
      'url': '/tests/joinedlist',
    },
  })
  static async TestsListJoined(): Promise<DG.DataFrame| undefined> {
    const pacakageTests = await TestAnalysisManager.collectPackageTests();
    const packageTestsListMapped = pacakageTests.map((elem) => {
      return {'type': 'package ', 'test': elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const manualTest = await TestAnalysisManager.collectManualTestNames();
    const manualTestsListMapped = manualTest.map((elem) => {
      return {'type': 'manual ', 'test': 'Test Track: ' + elem};
    });
    const resultTestsList = DG.DataFrame.fromObjects(manualTestsListMapped.concat(packageTestsListMapped));

    // const builds: DG.DataFrame = await grok.functions.call('UsageAnalysis:Builds');
    // const id = builds.get('name', 0);

    // const tests = await grok.functions.call('UsageAnalysis:getTestStatusesAcordingDF', { 'buildId': id, 'testslist': resultTestsList });
    grok.shell.addTableView(resultTestsList!);
    return resultTestsList;
  }

  @grok.decorators.func()
  static async TestAnalysisReportForCurrentDay(
    @grok.decorators.param({'type': 'datetime'}) date: any) : Promise<DG.DataFrame> {
    const tests = await TestAnalysisManager.collectPackageTests();
    const testsListMapped = tests.map((elem) => {
      return {'name': 'test-package ' + elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const testRuns = await grok.functions.call('UsageAnalysis:getServerStartTestResults', {'date': getDate(new Date(date)), 'testslist': DG.DataFrame.fromObjects(testsListMapped)});
    return testRuns;
  }

  @grok.decorators.app({
    'url': '/',
    'name': 'Usage Analysis',
    'meta': {'role': 'adminApp'},
  })
  static usageAnalysisApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string,
    @grok.decorators.param({'options': {'optional': true}}) date?: string,
    @grok.decorators.param({'options': {'optional': true}}) groups?: string,
    @grok.decorators.param({'options': {'optional': true}}) packages?: string,
    @grok.decorators.param({'options': {'optional': true}}) tags?: string,
    @grok.decorators.param({'options': {'optional': true}}) categories?: string,
    @grok.decorators.param({'options': {'optional': true}}) projects?: string): DG.ViewBase | null {
    const handler = new ViewHandler();
    handler.view.parentCall = grok.functions.getCurrentCall();
    handler.init(date, groups, packages, tags, categories, projects, path);
    return handler.view;
  }

  // Release-readiness dashboard. Dev-only: available and visible on dev.datagrok.ai; a no-op elsewhere.
  @grok.decorators.app({
    'url': '/release',
    'name': 'Release',
    'meta': {'role': 'adminApp'},
  })
  static releaseDashboardApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string): DG.ViewBase | null {
    if (!isDevHost()) {
      grok.shell.warning('The Release dashboard is available on dev.datagrok.ai only.');
      return null;
    }
    const handler = new ReleaseHandler(path);
    handler.view.parentCall = grok.functions.getCurrentCall();
    return handler.view;
  }

  /** Published VEX scan results (one row per scanned image) — feeds the
   * CicdTests dashboard's data-synced VEX table. Core images use their
   * bleeding-edge scan; packages/tools use the latest released scan. */
  @grok.decorators.func()
  static async vexImages(): Promise<DG.DataFrame> {
    return vexImagesToDataFrame(await fetchVexDashboardImages());
  }

  /** JS wrappers for the CicdTests dashboard's data-synced stress tables:
   * project datasync of DataQuery calls intermittently loses the result
   * dataframe (GROK-20391), while JS-function-backed sync is reliable. */
  @grok.decorators.func()
  static async stressSummaryDashboard(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:StressTestsSummary', {lastBuildsNum: 20});
  }

  @grok.decorators.func()
  static async stressRawDashboard(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:StressTestsRaw', {build: null});
  }

  @grok.decorators.app({
    'url': '/tests/manager',
    'name': 'Test Track',
    'icon': 'images/icons/test-track.svg',
    'meta': {'role': 'adminApp'},
  })
  static testTrackApp(): void {
    if (!grok.shell.dockManager.findNode(TestTrack.getInstance().root))
      TestTrack.getInstance().init();
    else
      TestTrack.getInstance().reopen();
  }


  @grok.decorators.app({
    'url': '/reports',
    'name': 'Reports',
    'icon': 'images/icons/reports.svg',
    'meta': {'role': 'adminApp', 'browsePath': 'Observability'},
  })
  static async reportsApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string): Promise<DG.ViewBase> {
    const parent = grok.functions.getCurrentCall();
    const app = new ReportingApp(parent);
    app.init(path).catch((e) => console.log(e));
    return app.view;
  }


  @grok.decorators.appTreeBrowser({app: 'Reports'})
  static async reportsAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    await treeNode.group('Reports', null, false).loadSources(grok.dapi.reports.by(10));
    await treeNode.group('Rules', null, false).loadSources(grok.dapi.rules.include('actions,actions.assignee').by(10));
  }

  @grok.decorators.dashboard({
    'meta': {
      'canView': 'Developers,Administrators',
    },
    'name': 'Usage',
    'test': 'usageWidget()',
  })
  static usageWidget(): DG.Widget {
    return new UsageWidget();
  }

  @grok.decorators.dashboard({
    'meta': {'canView': 'Developers,Administrators'},
    'name': 'Reports',
    'test': 'reportsWidget()',
  })
  static reportsWidget(): DG.Widget {
    return new ReportsWidget();
  }


  @grok.decorators.func({})
  static packageUsageWidget(
    @grok.decorators.param({'type': 'object'}) pckg: DG.Package): DG.Widget {
    return new PackageUsageWidget(pckg);
  }


  @grok.decorators.func({
    'meta': {showInGallery: 'false', role: 'viewer'},
    'outputs': [
      {
        type: 'viewer',
        name: 'result',
      },
    ],
  })
  static testDashboardsViewer(): TestDashboardWidget {
    return new TestDashboardWidget();
  }

  @grok.decorators.func({
    meta: {inspectorPanel: 'true'},
    name: 'Click Events',
  })
  static _specificClicksPanel(): DG.Widget {
    return new ClickEventsWidget();
  }

  @grok.decorators.autostart()
  static describeCurrentObj(): void {
    grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      const ent = acc.context;
      if (ent != null && ent.constructor.name === 'Package') {
        const pane = acc.addPane('Usage', () => ui.wait(async () => {
          let widget: HTMLElement;
          try {
            widget = PackageFunctions.packageUsageWidget(ent).root;
          } catch (e) {
            widget = ui.divText('Error on loading', {style: {color: 'var(--failure)'}});
          }
          return widget;
        }));
        const UAlink = ui.link('', async () => {
          grok.shell.v.path = `/apps/UsageAnalysis/Packages?date=this%20week&groups=${(await grok.dapi.groups.getGroupsLookup('All users'))[0].id}&packages=${ent.name}`;
          grok.functions.eval('UsageAnalysis:usageAnalysisApp()');
        }, 'Open Usage Analysis');
        UAlink.style.marginLeft = '3px';
        const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement;
        header.appendChild(UAlink);
      }
    });
  }
}
