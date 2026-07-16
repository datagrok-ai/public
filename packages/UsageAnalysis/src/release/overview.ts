import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {queries} from '../package-api';
import {fetchVexIndex, criticalHighCount, toDashboardImages, vexReleaseDelta} from '../tabs/vulnerabilities';
import {fetchReleaseTests, computeTestAlerts, stressRegression, defaultNextVersion, ReleaseContext, STALE_DAYS} from './data';
import {fetchReleaseTickets} from './tickets';

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

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Overview';
    this.ctx.env.subscribe(() => { if (this.initialized) this.refresh(); });
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.asOfLabel = ui.divText('Release readiness', 'ua-metrics-asof');
    this.bannerHost = ui.div([], 'ua-metrics-panel ua-release-alerts');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.asOfLabel, ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    const cardsRow = ui.divH([
      this.makeCard('Tests', 'pass rate', 'Tests'),
      this.makeCard('Failing', 'unmuted', 'Tests'),
      this.makeCard('Unstable', 'flaky tests', 'Tests'),
      this.makeCard('Stress', 'latest build', 'Stress'),
      this.makeCard('Vulnerabilities', 'critical + high', 'Vulnerabilities'),
      this.makeCard('Tickets', 'unresolved', 'Tickets'),
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

  private async refresh(): Promise<void> {
    this.bannerHost.innerHTML = '';
    this.bannerHost.append(ui.loader());
    const alerts: AlertNode[] = [];

    const [testsR, stressR, vexR, ticketsR] = await Promise.allSettled([
      fetchReleaseTests(this.ctx.env.value, 5),
      queries.stressTestsSummary(20),
      fetchVexIndex(),
      fetchReleaseTickets(defaultNextVersion()),
    ]);

    if (testsR.status === 'fulfilled' && testsR.value) {
      const a = computeTestAlerts(testsR.value);
      this.setCard('Tests', `${a.passRate}%`, a.passRate >= 95 ? 'green' : a.passRate >= 80 ? 'orange' : 'red',
        `${a.passed} passed, ${a.failed} failing`);
      this.setCard('Failing', `${a.failed}`, a.failed === 0 ? 'green' : 'red');
      this.setCard('Unstable', `${a.flaky}`, a.flaky === 0 ? 'green' : 'orange');
      if (a.failed > 0)
        alerts.push({label: `${a.failed} failing tests`, band: 'red', tab: 'Tests', groups: a.failingByPkg});
      if (a.slower > 0)
        alerts.push({label: `${a.slower} tests slower than recent median (passed runs)`, band: 'orange',
          tab: 'Tests', groups: a.slowerByPkg});
      if (a.flaky > 0)
        alerts.push({label: `${a.flaky} unstable/flaky tests`, band: 'orange', tab: 'Tests', groups: a.flakyByPkg});
      if (a.stale > 0)
        alerts.push({label: `${a.stale} tests not run for ${STALE_DAYS}+ days`, band: 'orange', tab: 'Tests', groups: a.staleByPkg});
    } else {
      this.setCard('Tests', 'n/a', 'info');
    }

    if (stressR.status === 'fulfilled') {
      const reg = stressRegression(stressR.value);
      this.setCard('Stress', reg.regressed ? 'Slower' : 'OK', reg.regressed ? 'red' : 'green');
      if (reg.regressed)
        alerts.push({label: `Stress slower — ${reg.detail}`, band: 'red', tab: 'Stress', groups: {}});
    } else {
      this.setCard('Stress', 'n/a', 'info');
    }

    if (vexR.status === 'fulfilled') {
      const {services, bleedingEdge} = vexR.value;
      const {critical, high} = criticalHighCount(toDashboardImages(services, bleedingEdge));
      this.setCard('Vulnerabilities', `${critical + high}`, critical > 0 ? 'red' : high > 0 ? 'orange' : 'green',
        `${critical} critical, ${high} high`);
      const delta = vexReleaseDelta(services, bleedingEdge);
      if (delta.newCritical > 0 || delta.newHigh > 0)
        alerts.push({label: `${delta.newCritical} new critical, ${delta.newHigh} new high vs last release`,
          band: delta.newCritical > 0 ? 'red' : 'orange', tab: 'Vulnerabilities',
          groups: {'': delta.images.map((x) => `${x.repo}: ${x.dCritical >= 0 ? '+' : ''}${x.dCritical} critical, ` +
            `${x.dHigh >= 0 ? '+' : ''}${x.dHigh} high`)}});
    } else {
      this.setCard('Vulnerabilities', 'n/a', 'info');
    }

    if (ticketsR.status === 'fulfilled') {
      const issues = ticketsR.value ?? [];
      const open = issues.filter((r: any) => (r.fields?.resolution?.name ?? 'Unresolved') === 'Unresolved').length;
      this.setCard('Tickets', `${open}`, open === 0 ? 'green' : 'info', `of ${issues.length} total`);
    } else {
      this.setCard('Tickets', 'n/a', 'info');
    }

    this.renderBanner(alerts);
  }

  private renderBanner(alerts: AlertNode[]): void {
    this.bannerHost.innerHTML = '';
    this.bannerHost.append(ui.divH([ui.divText('Alerts', 'ua-metrics-panel-title')], 'ua-metrics-panel-header'));
    if (alerts.length === 0) {
      this.bannerHost.append(ui.div([ui.divText('Release looks healthy — no blocking alerts.',
        'ua-metrics-card-sub ua-metrics-green')], 'ua-metrics-grid-host'));
      return;
    }
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
