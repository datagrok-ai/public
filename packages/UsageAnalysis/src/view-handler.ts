import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './tabs/ua';
import {UaToolbox} from './ua-toolbox';
import {EventsView} from './tabs/events';
import {PackagesView} from './tabs/packages';
import {FunctionsView} from './tabs/functions';
import {OverviewView} from './tabs/overview';
import {LogView} from './tabs/log';
import {SystemActivityView} from './tabs/system-activity';
import {ErrorsView} from './tabs/errors';
import {ProjectsView} from "./tabs/projects";
import {ClicksView} from './tabs/clicks';
import {MetricsView} from './tabs/metrics';
import {VulnerabilitiesView} from './tabs/vulnerabilities';
import {StressView} from './tabs/stress-tests';

export class ViewHandler {
  public static UA_NAME = 'Usage Analysis';
  private urlParams: Map<string, string> = new Map<string, string>();
  public view: DG.MultiView;

  constructor() {
    this.view = new DG.MultiView({viewFactories: {}});
    this.view.name = ViewHandler.UA_NAME;
    this.view.box = true;
  }

  async init(date?: string, groups?: string, packages?: string, tags?: string, categories?: string, projects?: string, path?: string): Promise<void> {
    const toolboxPromise = UaToolbox.construct(this);
    const viewClasses: (typeof UaView)[] = [OverviewView, PackagesView, FunctionsView, EventsView, ClicksView, LogView,
      SystemActivityView, ErrorsView, ProjectsView, MetricsView, StressView, VulnerabilitiesView];
    const views: UaView[] = [];
    for (let i = 0; i < viewClasses.length; i++) {
      const currentView = new viewClasses[i]();
      views.push(currentView);
      this.view.addView(currentView.name, () => {
        currentView.tryToInitViewers(path);
        return currentView;
      }, false);
    }

    let urlTab = 'Overview';

    if (path != undefined && path.length > 1) {
      const segments = path.split('/').filter((s) => s != '');
      const urlView = segments.length > 0 ? views.find((v) => ViewHandler.urlName(v) === segments[0].toLowerCase()) : undefined;
      if (urlView)
        urlTab = urlView.name;
    }

    const indicatorTimer = setTimeout(() => ui.setUpdateIndicator(this.view.root, true, 'Loading...'), 200);
    const toolbox = await toolboxPromise;
    clearTimeout(indicatorTimer);
    ui.setUpdateIndicator(this.view.root, false);
    for (const view of views)
      view.setToolbox(toolbox);
    this.view.toolbox = toolbox.rootAccordion.root;

    toolbox.toggleCategoriesInput(urlTab == 'Packages');
    toolbox.toggleTagsInput(urlTab == 'Functions');
    toolbox.toggleProjectsInput(urlTab == 'Projects');
    toolbox.togglePackagesInput(!['Projects', 'Metrics', 'Stress', 'Vulnerabilities', 'System Activity'].includes(urlTab));
    toolbox.toggleGroupsInput(!['Metrics', 'Stress', 'Vulnerabilities'].includes(urlTab));

    const paramsHaveDate = date != undefined;
    const paramsHaveUsers = groups != undefined;
    const paramsHavePackages = packages != undefined;
    const paramsHaveTags = tags != undefined;
    const paramsHavePackagesCategories = categories != undefined;
    const paramsHaveProjects = projects != undefined;
    if (paramsHaveDate || paramsHaveUsers || paramsHavePackages || paramsHavePackagesCategories || paramsHaveProjects) {
      if (paramsHaveDate)
        toolbox.setDate(date!);
      if (paramsHaveUsers)
        toolbox.setGroups(groups!);
      if (paramsHavePackages)
        toolbox.setPackages(packages!);
      if (paramsHaveTags)
        toolbox.setTags(tags!);
      if (paramsHavePackagesCategories)
        toolbox.setPackagesCategories(categories!);
      if (paramsHaveProjects)
        toolbox.setProjects(projects!);
      toolbox.applyFilter();
    }
    let helpShown = false;
    const puButton = ui.bigButton('Usage', () => {
      const v = this.getCurrentView();
      v.switchRout();
      this.updatePath();
      v.viewers[1].root.style.display = 'none';
      v.viewers[0].root.style.display = 'flex';
      puButton.disabled = true;
      piButton.disabled = false;
    });
    const piButton = ui.bigButton('Installation time', () => {
      const v = this.getCurrentView();
      v.switchRout();
      this.updatePath();
      v.viewers[0].root.style.display = 'none';
      v.viewers[1].root.style.display = 'flex';
      puButton.disabled = false;
      piButton.disabled = true;
    });
    puButton.disabled = true;
    const pButtons = ui.divH([puButton, piButton], 'ua-packages-buttons');
    pButtons.style.display = 'none';
    toolbox.filters.root.before(pButtons);

    const fuButton = ui.bigButton('Usage', () => {
      const v = this.getCurrentView() as FunctionsView;
      v.switchRout();
      this.updatePath();
      v.functionsExecTime.style.display = 'none';
      v.viewers[0].root.style.display = 'flex';
      fuButton.disabled = true;
      feButton.disabled = false;
    });
    const feButton = ui.bigButton('Execution time', () => {
      const v = this.getCurrentView() as FunctionsView;
      v.switchRout();
      this.updatePath();
      v.viewers[0].root.style.display = 'none';
      v.functionsExecTime.style.display = 'flex';
      fuButton.disabled = false;
      feButton.disabled = true;
    });
    fuButton.disabled = true;
    const fButtons = ui.divH([fuButton, feButton], 'ua-packages-buttons');
    fButtons.style.display = 'none';
    toolbox.filters.root.before(fButtons);

    this.view.tabs.onTabChanged.subscribe((_) => {
      const view = this.view.currentView;
      toolbox.toggleCategoriesInput(view.name === 'Packages');
      toolbox.toggleTagsInput(view.name === 'Functions');
      toolbox.toggleProjectsInput(view.name == 'Projects');
      toolbox.togglePackagesInput(!['Projects', 'Metrics', 'Stress', 'Vulnerabilities', 'System Activity'].includes(view.name));
      toolbox.toggleGroupsInput(!['Metrics', 'Stress', 'Vulnerabilities'].includes(view.name));
      // ViewHandler.UA.path = ViewHandler.UA.path.replace(/(UsageAnalysis\/)([a-zA-Z/]+)/, '$1' + view.name);
      this.updatePath();
      if (view instanceof UaView) {
        for (const viewer of view.viewers) {
          if (!viewer.activated) {
            viewer.activated = true;
            viewer.reloadViewer();
          }
        }
      }
      if (!helpShown) {
        if (this.view.currentView instanceof PackagesView || this.view.currentView instanceof FunctionsView) {
          grok.shell.windows.showToolbox = true;
          grok.shell.windows.showContextPanel = true;
          const info = ui.divText(`To learn more about an event, click the corresponding point.\
          To select multiple points, use CTRL + Click or SHIFT + Mouse Drag. Once you've made your selection,\
          see the detailed information on the Context Panel`);
          info.classList.add('ua-hint');
          grok.shell.o = info;
        }
        helpShown = true;
      }

      if (view.name === 'Packages')
        pButtons.style.display = 'flex';
      else
        pButtons.style.display = 'none';
      if (view.name === 'Functions')
        fButtons.style.display = 'flex';
      else
        fButtons.style.display = 'none';
    });
    this.changeTab(urlTab);
  }

  static urlName(view: DG.ViewBase): string {
    return view.name.replace(/\s/g, '').toLowerCase();
  }

  public getView(name: string) {
    return this.view.getView(name) as UaView;
  }

  public getCurrentView(): UaView {
    return this.view.currentView as UaView;
  }

  public changeTab(name: string) {
    this.view.tabs.currentPane = this.view.tabs.getPane(name);
  }

  setUrlParam(key: string, value: string, saveDuringChangingView: boolean = false) {
    const urlParams = new URLSearchParams(window.location.search);
    urlParams.set(key, value);

    const params: string[] = [];

    for (const keyAndValue of urlParams.entries())
      params.push(encodeURI(keyAndValue.join('=')));

    if (saveDuringChangingView)
      this.urlParams.set(key, value);

    this.view.path = `/${ViewHandler.urlName(this.getCurrentView())}?${params.join('&')}`.toLowerCase();
  }

   updatePath(): void {
    const v = this.getCurrentView();
    const s = this.view.path.split('?');
    const params = s.length === 2 ? s[1] : null;
     this.view.path = `/${ViewHandler.urlName(v)}${v.rout ?? ''}${params ? '?' + params : ''}`.toLowerCase();
  }
}
