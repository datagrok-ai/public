import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {fetchVexIndex, criticalHighCount, toDashboardImages, vexReleaseDelta} from '../tabs/vulnerabilities';
import {computeTestAlerts, fetchReleaseBenchmarks, computeBenchmarkAlerts, stressRegression, fetchStressSummary,
  ReleaseContext, ReleasePivot, BenchmarkPivot, STALE_DAYS, BENCH_SLOW_PCT} from './data';

type Band = 'green' | 'orange' | 'red' | 'info';
const BAND_CLASS: Record<Band, string> = {
  green: 'ua-metrics-green', orange: 'ua-metrics-orange', red: 'ua-metrics-red', info: 'ua-metrics-info'};

interface Card { root: HTMLElement; value: HTMLElement; sub: HTMLElement; }
interface AlertNode { label: string; band: Band; tab: string; groups: Record<string, string[]>; }

export class ReleaseOverviewView extends UaView {
  public switchTab: (name: string) => void = () => {};
  private asOfLabel!: HTMLElement;
  private cards = new Map<string, Card>();
  private bannerHost!: HTMLElement;
  private run = 0;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Overview';
    this.ctx.scopeChanged.subscribe(() => { if (this.initialized) this.refresh(); });
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.asOfLabel = ui.divText('Release readiness', 'ua-metrics-asof');
    this.bannerHost = ui.div([], 'ua-metrics-panel ua-release-alerts');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh(true));
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.asOfLabel, ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    const cardsRow = ui.divH([
      this.makeCard('Tests', 'pass rate', 'Tests'),
      this.makeCard('Failing', 'unmuted', 'Tests'),
      this.makeCard('Unstable', 'flaky tests', 'Tests'),
      this.makeCard('Benchmarks', 'slower than baseline', 'Benchmarks'),
      this.makeCard('Stress', 'latest build', 'Stress'),
      this.makeCard('Vulnerabilities', 'platform images', 'Vulnerabilities'),
    ], 'ua-metrics-cards-row');

    this.root.append(ui.div([header, cardsRow, this.bannerHost], 'ua-metrics-root'));
    await this.refresh();
  }

  private makeCard(title: string, sub: string, target: string): HTMLElement {
    const value = ui.divText('—', 'ua-metrics-card-value');
    const subEl = ui.divText(sub, 'ua-metrics-card-sub');
    const root = ui.div([ui.divText(title, 'ua-metrics-card-title'), value, subEl], 'ua-metrics-card');
    root.style.cursor = 'pointer';
    root.onclick = () => this.switchTab(target);
    this.cards.set(title, {root, value, sub: subEl});
    return root;
  }

  private setCard(title: string, value: string, band: Band, sub?: string): void {
    const c = this.cards.get(title);
    if (!c) return;
    c.value.textContent = value;
    c.value.className = `ua-metrics-card-value ${BAND_CLASS[band]}`;
    if (sub != null)
      c.sub.textContent = sub;
  }

  /** Every source loads and renders on its own: a slow or failing one holds back only its card. */
  private async refresh(force = false): Promise<void> {
    const run = ++this.run;
    const sections: AlertNode[][] = [];
    let pending = 0;
    const load = async <T>(cards: string[], data: Promise<T>, apply: (v: T) => AlertNode[]): Promise<void> => {
      const i = pending++;
      let alerts: AlertNode[] = [];
      try {
        const v = await data;
        if (run === this.run)
          alerts = apply(v);
      } catch (e) {
        grok.log.warning(`Release overview: ${cards[0]} failed: ${e}`);
        if (run === this.run)
          cards.forEach((c) => this.setCard(c, 'n/a', 'info'));
      }
      if (run !== this.run)
        return;
      sections[i] = alerts;
      this.renderBanner(sections.flat(), --pending > 0);
    };
    const all = [
      load(['Tests', 'Failing', 'Unstable'], this.ctx.tests(5, force), (p) => this.applyTests(p)),
      load(['Benchmarks'], fetchReleaseBenchmarks(this.ctx.env.value, 5, this.ctx.branch.value),
        (p) => this.applyBenchmarks(p)),
      load(['Stress'], fetchStressSummary(20), (df) => this.applyStress(df)),
      load(['Vulnerabilities'], fetchVexIndex(), (vex) => this.applyVex(vex)),
    ];
    this.renderBanner([], true);
    await Promise.all(all);
  }

  private applyTests(p: ReleasePivot | null): AlertNode[] {
    if (!p) {
      ['Tests', 'Failing', 'Unstable'].forEach((c) => this.setCard(c, 'n/a', 'info'));
      return [];
    }
    const alerts: AlertNode[] = [];
    const a = computeTestAlerts(p);
    this.setCard('Tests', `${a.passRate}%`, a.passRate >= 95 ? 'green' : a.passRate >= 80 ? 'orange' : 'red',
      `${a.passed} passed, ${a.failed} failing`);
    this.setCard('Failing', `${a.failed}`, a.failed === 0 ? 'green' : 'red');
    this.setCard('Unstable', `${a.flaky}`, a.flaky === 0 ? 'green' : 'orange');
    if (a.failed > 0)
      alerts.push({label: `${a.failed} failing tests`, band: 'red', tab: 'Tests', groups: a.failingByPkg});
    if (a.slower > 0)
      alerts.push({label: `${a.slower} tests slower than baseline (recent avg or previous release)`,
        band: 'orange', tab: 'Tests', groups: a.slowerByPkg});
    if (a.flaky > 0)
      alerts.push({label: `${a.flaky} unstable/flaky tests`, band: 'orange', tab: 'Tests', groups: a.flakyByPkg});
    if (a.stale > 0)
      alerts.push({label: `${a.stale} tests not run for ${STALE_DAYS}+ days`, band: 'orange', tab: 'Tests', groups: a.staleByPkg});
    return alerts;
  }

  // A slower benchmark is not release-blocking on its own (orange, not red), but a benchmark that
  // outright FAILED blew its own threshold inside the suite — that is a red.
  private applyBenchmarks(p: BenchmarkPivot | null): AlertNode[] {
    if (!p) {
      this.setCard('Benchmarks', 'n/a', 'info', `no runs within ${BENCH_SLOW_PCT}% baseline yet`);
      return [];
    }
    const alerts: AlertNode[] = [];
    const b = computeBenchmarkAlerts(p);
    this.setCard('Benchmarks', b.slower === 0 ? 'OK' : `${b.slower}`,
      b.failed > 0 ? 'red' : b.slower === 0 ? 'green' : 'orange',
      b.slower === 0 ? `${b.total} within ${p.slowPct}%` : `worst +${b.worstPct}% of ${b.total}`);
    if (b.slower > 0)
      alerts.push({label: `${b.slower} benchmarks slower than expected (>${p.slowPct}% over baseline)`,
        band: 'orange', tab: 'Benchmarks', groups: b.slowerByPkg});
    if (b.failed > 0)
      alerts.push({label: `${b.failed} benchmarks failed`, band: 'red', tab: 'Benchmarks', groups: b.failingByPkg});
    return alerts;
  }

  private applyStress(df: DG.DataFrame): AlertNode[] {
    const reg = stressRegression(df);
    this.setCard('Stress', reg.regressed ? 'Slower' : 'OK', reg.regressed ? 'red' : 'green');
    return reg.regressed ? [{label: `Stress slower — ${reg.detail}`, band: 'red', tab: 'Stress', groups: {}}] : [];
  }

  private applyVex({services, bleedingEdge}: Awaited<ReturnType<typeof fetchVexIndex>>): AlertNode[] {
    // The card answers for the platform's own images (what the next release ships). Plugin and
    // tool images release on their own schedule; they are counted in the sub-line and on the tab.
    const {critical, high} = criticalHighCount(bleedingEdge.filter((r) => r.category === 'core'));
    const all = criticalHighCount(toDashboardImages(services, bleedingEdge));
    this.setCard('Vulnerabilities', `${critical + high}`, critical > 0 ? 'red' : high > 0 ? 'orange' : 'green',
      `${critical} critical, ${high} high; ${all.critical + all.high - critical - high} in plugins`);
    const delta = vexReleaseDelta(services, bleedingEdge);
    if (delta.newCritical === 0 && delta.newHigh === 0)
      return [];
    return [{label: `${delta.newCritical} new critical, ${delta.newHigh} new high vs last release`,
      band: delta.newCritical > 0 ? 'red' : 'orange', tab: 'Vulnerabilities',
      groups: {'': delta.images.map((x) => `${x.repo}: ${x.dCritical >= 0 ? '+' : ''}${x.dCritical} critical, ` +
        `${x.dHigh >= 0 ? '+' : ''}${x.dHigh} high`)}}];
  }

  private renderBanner(alerts: AlertNode[], loading: boolean): void {
    this.bannerHost.innerHTML = '';
    this.bannerHost.append(ui.divH([ui.divText('Alerts', 'ua-metrics-panel-title')], 'ua-metrics-panel-header'));
    if (loading)
      this.bannerHost.append(ui.loader());
    else if (alerts.length === 0) {
      this.bannerHost.append(ui.div([ui.divText('Release looks healthy — no blocking alerts.',
        'ua-metrics-card-sub ua-metrics-green')], 'ua-metrics-grid-host'));
    }
    if (alerts.length === 0)
      return;
    const tree = ui.tree();
    tree.root.classList.add('ua-release-alert-tree');
    for (const a of alerts) {
      const total = Object.values(a.groups).reduce((n, xs) => n + xs.length, 0);
      if (total === 0) {
        const leaf = tree.item(a.label);
        leaf.captionLabel.classList.add(BAND_CLASS[a.band]);
        leaf.root.style.cursor = 'pointer';
        leaf.root.onclick = () => this.switchTab(a.tab);
        continue;
      }
      const cat = tree.group(`${a.label} (${total})`, null, false);
      cat.captionLabel.classList.add(BAND_CLASS[a.band]);
      // Second level: group the tests by package (an empty package key lists items directly).
      for (const pkg of Object.keys(a.groups).sort()) {
        const items = a.groups[pkg];
        if (!items.length)
          continue;
        const parent = pkg === '' ? cat : cat.group(`${pkg} (${items.length})`, null, false);
        for (const it of items.slice(0, 200)) {
          const node = parent.item(it);
          node.root.style.cursor = 'pointer';
          node.root.onclick = () => this.switchTab(a.tab);
        }
        if (items.length > 200)
          parent.item(`… and ${items.length - 200} more`);
      }
    }
    this.bannerHost.append(ui.div([tree.root], 'ua-metrics-grid-host'));
  }
}
