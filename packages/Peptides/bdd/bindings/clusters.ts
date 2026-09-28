import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {el, expect, locate, viewers, type ElementRef} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

type ClusterStats = {count: number; mean: number; meanDifference: number; hasPValue: boolean};
type ClusterSource = {total: number; clusters: Record<string, ClusterStats>};

/** Independent arithmetic over the analysis table, without the model's statistics or selection. */
async function sourceClusters(page: Page, source: string): Promise<ClusterSource> {
  return page.evaluate((name) => {
    const table = grok.shell.tables.find((t: any) => t.name === name);
    if (!table || table.rowCount === 0)
      throw new Error(`no populated source table "${name}"`);
    const clusterColumn = table.getCol('Cluster (MCL)');
    const activityColumn = table.getCol('Activity');
    const groups = new Map<string, number[]>();
    const activities: {cluster: string; activity: number}[] = [];
    for (let row = 0; row < table.rowCount; row++) {
      if (!table.filter.get(row))
        continue;
      const cluster = String(clusterColumn.get(row));
      const activity = activityColumn.get(row);
      if (clusterColumn.isNone(row) || activityColumn.isNone(row) || !Number.isFinite(activity))
        throw new Error(`source row ${row + 1} lacks a cluster or finite Activity`);
      activities.push({cluster, activity});
      const group = groups.get(cluster) ?? [];
      group.push(activity);
      groups.set(cluster, group);
    }
    if (!activities.length)
      throw new Error(`source table "${name}" has no filtered rows`);
    const clusters: Record<string, ClusterStats> = {};
    for (const [cluster, selected] of groups) {
      const rest = activities.filter((row) => row.cluster !== cluster).map((row) => row.activity);
      const mean = selected.reduce((sum, value) => sum + value, 0) / selected.length;
      const restMean = rest.reduce((sum, value) => sum + value, 0) / Math.max(rest.length, 1);
      clusters[cluster] = {count: selected.length, mean, meanDifference: mean - restMean,
        hasPValue: selected.length > 1 && rest.length > 1};
    }
    return {total: activities.length, clusters};
  }, source);
}

export const clusterSummaries = Then('the cluster summaries should match the source rows of table {string}',
  async (page: Page, source: string) => {
    const expected = await sourceClusters(page, source);
    const target = el('Logo Summary Table viewer');
    await viewers.settle(page, target);
    const values = await viewers.onViewer(page, target, (root) =>
      (window as any).__bdd.viewerOf(root).getWidgetStatus().values);
    expect(values.clusters, 'one summary row for every source cluster').toBe(Object.keys(expected.clusters).length);
    expect(values['members total'], 'cluster members account for every source row').toBe(expected.total);
    for (const [cluster, stats] of Object.entries(expected.clusters)) {
      expect(values[`Members of cluster ${cluster}`], `source members of cluster ${cluster}`).toBe(stats.count);
      const difference = values[`Mean difference of cluster ${cluster}`];
      expect(Number.isFinite(difference), `finite mean difference of cluster ${cluster}`).toBe(true);
      // The displayed summary columns are Float32; arithmetic above uses the source Activity values.
      expect(Math.abs(difference - stats.meanDifference), `mean difference of cluster ${cluster} against the rest`)
        .toBeLessThanOrEqual(Math.max(Math.abs(stats.meanDifference) * 1e-6, 1e-10));
      const pValue = values[`P-Value of cluster ${cluster}`];
      if (stats.hasPValue) {
        expect(Number.isFinite(pValue), `finite p-value of cluster ${cluster}`).toBe(true);
        expect(pValue, `p-value of cluster ${cluster}`).toBeGreaterThanOrEqual(0);
        expect(pValue, `p-value of cluster ${cluster}`).toBeLessThanOrEqual(1);
      } else
        expect(pValue, `no t-test when cluster ${cluster} or its complement has fewer than two rows`).toBe('');
    }
  }, {description: 'all LST membership counts and mean differences against source Cluster (MCL)/Activity rows; valid p-values where a t-test is defined'});

export const clusterStatistics = Then('the cluster statistics in {element} should match cluster {string} of table {string}',
  async (page: Page, target: ElementRef, cluster: string, source: string) => {
    const expected = await sourceClusters(page, source);
    const stats = expected.clusters[cluster];
    if (!stats)
      throw new Error(`source table "${source}" has no cluster "${cluster}"; it has: ${Object.keys(expected.clusters).join(', ')}`);
    const row = async (label: string) =>
      (await locate(page, `"${label}" table row in ${target.phrase}`)).filter({visible: true});
    const displayed = {
      Count: `${stats.count} (${(stats.count / expected.total * 100).toFixed(3)}%)`,
      'Mean difference': stats.meanDifference.toFixed(3),
      'Mean activity': stats.mean.toFixed(3),
    };
    for (const [label, value] of Object.entries(displayed)) {
      const statRow = await row(label);
      await expect(statRow, `${label} row in ${target.phrase}`).toHaveCount(1);
      await expect(statRow.locator('td').nth(1), `${label} for source cluster ${cluster}`).toHaveText(value);
    }
    const pValueRow = await row('p-value');
    await expect(pValueRow, `p-value availability for source cluster ${cluster}`).toHaveCount(stats.hasPValue ? 1 : 0);
    if (stats.hasPValue)
      await expect(pValueRow.locator('td').nth(1), `p-value for source cluster ${cluster} is in [0, 1]`)
        .toHaveText(/^(?:<0\.01|0(?:\.\d+)?|1(?:\.0+)?)$/);
  }, {description: 'visible Count, Mean difference and Mean activity from independent source-row arithmetic; p-value range only, without assuming stable MCL membership'});
