import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {StressView} from '../tabs/stress-tests';
import {VulnerabilitiesView} from '../tabs/vulnerabilities';
import {ReleaseOverviewView} from './overview';
import {TestsView} from './tests';
import {ManualTestsView} from './manual';
import {ReleaseTicketsView} from './tickets';
import {ReleaseContext, ENV_CHOICES, DEFAULT_ENV} from './data';

// Assembles the focused /release dashboard: a MultiView of Overview + Tests + Stress +
// Vulnerabilities + Tickets. Stress and Vulnerabilities are reused verbatim from the main app; all
// tabs are toolbox-independent, so we skip the shared UaToolbox (markToolboxReady unblocks init).
export class ReleaseHandler {
  public static NAME = 'Release';
  public view: DG.MultiView;
  private ctx = new ReleaseContext();
  private envInput: DG.InputBase;
  private refreshIcon: HTMLElement;

  constructor(path?: string) {
    this.view = new DG.MultiView({viewFactories: {}});
    this.view.name = ReleaseHandler.NAME;
    this.view.box = true;

    this.envInput = ui.input.choice('Environment', {value: this.ctx.env.value, items: ENV_CHOICES,
      onValueChanged: () => this.ctx.env.next(this.envInput.value ?? DEFAULT_ENV)});
    this.envInput.setTooltip('Filter builds by the instance they ran on (applies to Overview and Tests)');
    this.refreshIcon = ui.iconFA('sync-alt', () => this.ctx.refresh.next(), 'Refresh the current data');
    this.refreshIcon.classList.add('ua-release-ribbon-refresh');

    const overview = new ReleaseOverviewView(this.ctx);
    overview.switchTab = (name) => this.changeTab(name);
    const views: UaView[] = [overview, new TestsView(this.ctx), new ManualTestsView(this.ctx),
      new StressView(undefined, this.ctx), new VulnerabilitiesView(undefined, this.ctx),
      new ReleaseTicketsView(this.ctx)];
    for (const v of views) {
      v.markToolboxReady();
      this.view.addView(v.name, () => {
        v.tryToInitViewers(path);
        return v;
      }, false);
    }

    // MultiView overwrites its ribbon with the active child's on every tab switch; child tabs keep their
    // controls in their own in-view header (empty ribbon), so re-assert the global env picker each time.
    this.view.tabs.onTabChanged.subscribe(() => {
      this.view.path = `/${this.view.currentView.name}`.toLowerCase();
      this.view.setRibbonPanels([[this.envInput.root, this.refreshIcon]]);
    });

    let startTab = 'Overview';
    if (path && path.length > 1) {
      const seg = path.split('/').filter((s) => s !== '')[0];
      if (seg)
        startTab = seg[0].toUpperCase() + seg.slice(1);
    }
    if (views.some((v) => v.name === startTab))
      this.changeTab(startTab);
    this.view.setRibbonPanels([[this.envInput.root, this.refreshIcon]]);
  }

  changeTab(name: string): void {
    this.view.tabs.currentPane = this.view.tabs.getPane(name);
  }
}
