import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {StressView} from '../tabs/stress-tests';
import {VulnerabilitiesView} from '../tabs/vulnerabilities';
import {ReleaseOverviewView} from './overview';
import {TestsView} from './tests';
import {BenchmarksView} from './benchmarks';
import {ManualTestsView} from './manual';
import {ReleaseTicketsView} from './tickets';
import {ReleaseContext, ENV_CHOICES, DEFAULT_ENV, DEFAULT_BRANCH, ALL_BRANCHES, fetchTestBranches,
  defaultBranchForEnv} from './data';

// Assembles the focused /release dashboard: a MultiView of Overview + Tests + Stress +
// Vulnerabilities + Tickets. Stress and Vulnerabilities are reused verbatim from the main app; all
// tabs are toolbox-independent, so we skip the shared UaToolbox (markToolboxReady unblocks init).
export class ReleaseHandler {
  public static NAME = 'Release';
  public view: DG.MultiView;
  private ctx = new ReleaseContext();
  private envInput: DG.InputBase;
  private branchInput: DG.InputBase;
  private branches: string[] = [];
  private refreshIcon: HTMLElement;

  constructor(path?: string) {
    this.view = new DG.MultiView({viewFactories: {}});
    this.view.name = ReleaseHandler.NAME;
    this.view.box = true;

    this.envInput = ui.input.choice('Environment', {value: this.ctx.env.value, items: ENV_CHOICES,
      onValueChanged: () => this.selectEnv(this.envInput.value ?? DEFAULT_ENV)});
    this.envInput.setTooltip('Filter builds by the instance they ran on (applies to Overview and Tests)');
    // Scoped to one branch by default: an ad-hoc run on a branch used to land in the same column
    // as the nightly and read as trunk (GROK-20760). Which branch follows the instance — dev runs
    // trunk, the others run a release — and switching instance re-picks it. The list is filled in
    // once the branches query answers; until then the current value is the only choice.
    this.branchInput = ui.input.choice('Branch', {value: this.ctx.branch.value,
      items: [this.ctx.branch.value, ALL_BRANCHES],
      onValueChanged: () => this.ctx.branch.next(this.branchInput.value ?? DEFAULT_BRANCH)});
    this.branchInput.setTooltip('Filter builds by the branch they were produced from (applies to Overview and Tests)');
    fetchTestBranches()
      .then((branches) => {
        this.branches = branches;
        (this.branchInput as DG.ChoiceInput<string>).items = branches;
        this.selectEnv(this.ctx.env.value);
      })
      .catch((e) => grok.log.warning(`Release: could not list test branches: ${e}`));
    this.refreshIcon = ui.iconFA('sync-alt', () => this.ctx.refresh.next(), 'Refresh the current data');
    this.refreshIcon.classList.add('ua-release-ribbon-refresh');

    const overview = new ReleaseOverviewView(this.ctx);
    overview.switchTab = (name) => this.changeTab(name);
    const views: UaView[] = [overview, new TestsView(this.ctx), new BenchmarksView(this.ctx),
      new ManualTestsView(this.ctx), new StressView(undefined, this.ctx),
      new VulnerabilitiesView(undefined, this.ctx), new ReleaseTicketsView(this.ctx)];
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
      this.view.setRibbonPanels([[this.envInput.root, this.branchInput.root, this.refreshIcon]]);
    });

    let startTab = 'Overview';
    if (path && path.length > 1) {
      const seg = path.split('/').filter((s) => s !== '')[0];
      if (seg)
        startTab = seg[0].toUpperCase() + seg.slice(1);
    }
    if (views.some((v) => v.name === startTab))
      this.changeTab(startTab);
    this.view.setRibbonPanels([[this.envInput.root, this.branchInput.root, this.refreshIcon]]);
  }

  /** Switches instance and re-points the branch at the one that instance runs. Both land in the
   * same tick, so the tabs coalesce them into a single re-fetch. */
  private selectEnv(env: string): void {
    const branch = defaultBranchForEnv(env, this.branches);
    if (this.branchInput.value !== branch)
      this.branchInput.value = branch;
    this.ctx.env.next(env);
    if (this.ctx.branch.value !== branch)
      this.ctx.branch.next(branch);
  }

  changeTab(name: string): void {
    this.view.tabs.currentPane = this.view.tabs.getPane(name);
  }
}
