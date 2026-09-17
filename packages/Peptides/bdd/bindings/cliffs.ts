import {expect, type Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {el, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

export const cliffParticipants = Then('the position {int} cliff chart should contain exactly its participating peptides',
  async (page: Page, position: number) => {
    const target = el('line chart viewer in Sequence Mutation Cliffs viewer');
    await viewers.settle(page, target);
    const data = await viewers.onViewer(page, target, (root, position) => {
      const source = grok.shell.t;
      const sequences = source.getCol('AlignedSequence');
      if (sequences.meta.units !== 'separator' || !sequences.getTag('separator'))
        throw new Error('the cliff-chart oracle requires separator notation');
      const chart = DG.Widget.find(root);
      return {
        split: sequences.toList().map((s: string) => s.split(sequences.getTag('separator'))) as string[][],
        activities: source.getCol('IC50').toList() as number[],
        // a missing monomer is a hole in the list, which map() would skip
        monomers: Array.from(chart.dataFrame.getCol(`Position ${position}`).toList(), (m) => (m as string) ?? ''),
        plotted: chart.dataFrame.getCol('IC50').toList() as number[],
      };
    }, position);
    expect(position).toBeGreaterThan(0);
    expect(data.split.every((s) => s.length === data.split[0].length && position <= s.length)).toBe(true);
    const participants = new Set<number>();
    for (let a = 0; a < data.split.length; a++) {
      for (let b = a + 1; b < data.split.length; b++) {
        if (data.split[a][position - 1] !== data.split[b][position - 1] &&
          data.split[a].every((m, p) => p === position - 1 || m === data.split[b][p])) {
          participants.add(a);
          participants.add(b);
        }
      }
    }
    expect(participants.size, 'the fixture has participating peptides').toBeGreaterThan(0);
    // The chart stores Float32 activities; these remain unique identifiers in this fixture.
    expect(new Set(data.activities.map(Math.fround)).size).toBe(data.activities.length);
    const expected = [...participants].sort((a, b) => a - b)
      .map((row) => [data.split[row][position - 1], Math.fround(data.activities[row])]);
    expect(data.monomers.map((monomer, row) => [monomer, Math.fround(data.plotted[row])]))
      .toEqual(expected);
  }, {description: 'compares chart monomers and activities with source pairs differing only at the requested position'});
