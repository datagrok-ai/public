/* eslint-disable max-len */
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {CatalogTable, Scope} from '@datagrok-libraries/u2';
import {CampaignGrouping, i18n} from '../app/consts';
import {AnyCampaign, CampaignFamily, catalogColumns, catalogGroupings, dayText, fieldText, fieldValue}
  from '../app/shared/campaign-catalog';
import {HitDesignCampaign} from '../app/types';
import {expectSame, formDate, oldDate} from './fixtures';

/** The fixed record list: every value shape the landing meets (codes, YYYY/M/D dates, counts, blanks, Date fields). */
function records(): HitDesignCampaign[] {
  const r = (x: Partial<HitDesignCampaign> & {name: string}) => ({templateName: 'Kinase', status: 'In Progress',
    createDate: '2024/1/1', campaignFields: {}, ...x}) as HitDesignCampaign;
  return [
    r({name: 'KHL-2', friendlyName: 'Kinase b', template: {key: 'KHL', name: 'Kinase'} as HitDesignCampaign['template'],
      createDate: '2024/10/17', rowCount: 12, authorUserFriendlyName: 'Ann', lastModifiedUserName: 'Bob',
      campaignFields: {'Chemist': 'Ann', 'Target MW': 350, 'Deadline': formDate}}),
    r({name: 'KHL-12', friendlyName: 'Kinase a', template: {key: 'KHL', name: 'Kinase'} as HitDesignCampaign['template'],
      status: 'Approved', createDate: '2024/9/4', rowCount: 3, authorUserFriendlyName: 'Cy', lastModifiedUserName: 'Ann',
      campaignFields: {'Chemist': 'Cy', 'Target MW': 90}}),
    r({name: 'GPCR-1', templateName: 'GPCR series', status: '', createDate: '2023/11/11'}),
    r({name: 'DMT1-1', friendlyName: 'Demo 1', template: {key: 'DMT1', name: 'Demo Template 1'} as HitDesignCampaign['template'],
      createDate: '2026/6/17', rowCount: 100, authorUserFriendlyName: 'Bob', lastModifiedUserName: 'Bob',
      campaignFields: {'Chemist': 'Ann', 'Target MW': '1000', 'Deadline': oldDate}}),
  ];
}

/** The order the campaigns table shows after a sort by the column, its own comparator. */
function sorted(columnName: string, ascending: boolean, list: AnyCampaign[] = records(), family: CampaignFamily = 'HD'): string[] {
  const column = catalogColumns(family, list, () => 0).find((c) => c.name === columnName)!;
  return CatalogTable.sortItems([...list], column, ascending).map((c) => c.name);
}

/** The table's groups for a grouping, by label in the order shown (a Map: an object would put "90" before "350"),
 * each with its rows in order. */
function groups(groupBy: string | null, list: AnyCampaign[] = records()): Map<string, string[]> {
  const scope = new Scope();
  try {
    const table = scope.build(() => new CatalogTable<AnyCampaign>({items: list, keyOf: (c) => c.name,
      columns: catalogColumns('HD', list, () => 0), groupings: catalogGroupings(list),
      emptyGroupLabel: i18n.noInformation, singleGroupHeader: true, state: {groupBy}}));
    const out = new Map<string, string[]>();
    let label = '';
    for (const item of table.getItems({all: true}).list) {
      if (item.key.startsWith('section:'))
        out.set(label = item.label, []);
      else
        out.set(label, [...out.get(label) ?? [], item.key]);
    }
    return out;
  } finally {
    scope.dispose();
  }
}

const members = (groupBy: string | null) => Object.fromEntries(groups(groupBy));
const labels = (groupBy: string) => [...groups(groupBy).keys()];

category('Catalog', () => {
  test('sort: created dates in YYYY/M/D sort as dates', async () => {
    expectSame(sorted('Created', true), ['GPCR-1', 'KHL-12', 'KHL-2', 'DMT1-1'], 'ascending');
    expectSame(sorted('Created', false), ['DMT1-1', 'KHL-2', 'KHL-12', 'GPCR-1'], 'descending');
  });

  test('sort: molecule counts sort as numbers, a missing count as 0', async () => {
    expectSame(sorted('Molecules', true), ['GPCR-1', 'KHL-12', 'KHL-2', 'DMT1-1'], 'ascending');
    expectSame(sorted('Molecules', false), ['DMT1-1', 'KHL-2', 'KHL-12', 'GPCR-1'], 'descending');
  });

  test('sort: codes, names, statuses and authors sort as text, blanks last both ways, ties stable', async () => {
    expectSame(sorted('Code', true), ['DMT1-1', 'GPCR-1', 'KHL-2', 'KHL-12'], 'Code ascending');
    expectSame(sorted('Name', true), ['DMT1-1', 'GPCR-1', 'KHL-12', 'KHL-2'], 'Name (friendly name, else code) ascending');
    expectSame(sorted('Status', true), ['KHL-12', 'KHL-2', 'DMT1-1', 'GPCR-1'], 'Status ascending');
    expectSame(sorted('Status', false), ['KHL-2', 'DMT1-1', 'KHL-12', 'GPCR-1'], 'Status descending');
    expectSame(sorted('Author', true), ['KHL-2', 'DMT1-1', 'KHL-12', 'GPCR-1'], 'Author ascending');
    expectSame(sorted('Last Modified by', false), ['KHL-2', 'DMT1-1', 'KHL-12', 'GPCR-1'], 'Last Modified by descending');
  });

  test('sort: codes sort in natural order and are never read as dates', async () => {
    const list = (names: string[]) => names.map((name) => ({name, campaignFields: {}}) as unknown as HitDesignCampaign);
    expectSame(sorted('Code', true, list(['KHL-2', 'KHL-12', 'KHL-1'])), ['KHL-1', 'KHL-2', 'KHL-12'], 'KHL-2 before KHL-12');
    expectSame(sorted('Code', true, list(['KHL-2', 'KHL-12', 'KHL-13', 'KHL-1'])), ['KHL-1', 'KHL-2', 'KHL-12', 'KHL-13'],
      'a code past 12 changes nothing');
  });

  test('sort: campaign fields, numbers among them, blanks last', async () => {
    expectSame(sorted('campaignFields.Target MW', true), ['KHL-12', 'KHL-2', 'DMT1-1', 'GPCR-1'], 'numbers ascending, blank last');
    expectSame(sorted('campaignFields.Target MW', false), ['DMT1-1', 'KHL-2', 'KHL-12', 'GPCR-1'], 'numbers descending, blank last');
    expectSame(sorted('campaignFields.Chemist', true), ['KHL-2', 'DMT1-1', 'KHL-12', 'GPCR-1'], 'text ascending, blank last');
  });

  test('sort: a Date campaign field orders by the date, in either stored shape', async () => {
    expectSame(sorted('campaignFields.Deadline', true), ['DMT1-1', 'KHL-2', 'KHL-12', 'GPCR-1'], 'DMT1-1 (2025-10-08) before KHL-2 (2026-03-27)');
  });

  test('sort: the Hit Triage columns sort Total and Selected', async () => {
    const triage = records().map((c, i) => ({...c, filteredRowCount: [5, 1, 0, 2][i]})) as unknown as AnyCampaign[];
    expectSame(sorted('Selected', false, triage, 'HT'), ['KHL-2', 'DMT1-1', 'KHL-12', 'GPCR-1'], 'Selected descending');
    expectSame(sorted('Total', true, triage, 'HT'), ['GPCR-1', 'KHL-12', 'KHL-2', 'DMT1-1'], 'Total ascending');
  });

  test('group: by template name, status, author, last modifier, and none; the empty group last', async () => {
    expectSame(members(null), {'': ['KHL-2', 'KHL-12', 'GPCR-1', 'DMT1-1']}, 'None');
    expectSame(labels(CampaignGrouping.Template), ['Kinase', 'GPCR series', 'Demo Template 1'],
      'Template: the template\'s name, in order of first appearance');
    expectSame(members(CampaignGrouping.Template), {'Kinase': ['KHL-2', 'KHL-12'], 'GPCR series': ['GPCR-1'], 'Demo Template 1': ['DMT1-1']}, 'Template members');
    expectSame(members(CampaignGrouping.Status), {'In Progress': ['KHL-2', 'DMT1-1'], 'Approved': ['KHL-12'], [i18n.noInformation]: ['GPCR-1']}, 'Status');
    expectSame(labels(CampaignGrouping.Author), ['Ann', 'Cy', 'Bob', i18n.noInformation], 'Author: first appearance, the empty last');
    expectSame(members(CampaignGrouping.LastModifiedUser), {'Bob': ['KHL-2', 'DMT1-1'], 'Ann': ['KHL-12'], [i18n.noInformation]: ['GPCR-1']},
      'Last Modified User');
  });

  test('group: a template is one group by its name, with or without a snapshot, whatever its spacing or case', async () => {
    const legacy = {name: 'KHL-30', templateName: ' kinase ', status: 'In Progress', createDate: '2024/1/1', campaignFields: {}} as
      unknown as HitDesignCampaign;
    expectSame(Object.fromEntries(groups(CampaignGrouping.Template, [...records(), legacy])),
      {'Kinase': ['KHL-2', 'KHL-12', 'KHL-30'], 'GPCR series': ['GPCR-1'], 'Demo Template 1': ['DMT1-1']},
      'the snapshot-less KHL-30 joins Kinase, under the first campaign\'s spelling');
  });

  test('group: the menu offers a grouping under the caption of the column it groups by', async () => {
    const shown = (name: string) => catalogGroupings(records()).find((g) => g.name === name)?.label ?? name;
    const columns = catalogColumns('HD', records(), () => 0).map((c) => c.header ?? c.name);
    for (const g of [CampaignGrouping.Status, CampaignGrouping.Author, CampaignGrouping.LastModifiedUser])
      expect(columns.includes(shown(g)), true, `"${shown(g)}" is a column's caption (${columns.join(', ')})`);
    expect(shown(CampaignGrouping.LastModifiedUser), 'Last Modified by', 'Last Modified User is offered as Last Modified by');
  });

  test('group: by a campaign field, in order of first appearance', async () => {
    expectSame(members('campaignFields.Chemist'), {'Ann': ['KHL-2', 'DMT1-1'], 'Cy': ['KHL-12'], [i18n.noInformation]: ['GPCR-1']}, 'Chemist');
    expectSame(labels('campaignFields.Target MW'), ['350', '90', '1000', i18n.noInformation], 'Target MW');
  });

  test('group: a Date campaign field gives one group per date, labelled as its cells show the date', async () => {
    const keys = labels('campaignFields.Deadline');
    expect(!keys.some((k) => k.includes('[object Object]')) && keys.length === 3, true,
      `groups for 2026-03-27, 2025-10-08 and no information; got ${JSON.stringify(keys)}`);
    expectSame(keys.filter((k) => k !== i18n.noInformation),
      [fieldText(formDate), fieldText(oldDate)].filter((d, i, all) => all.indexOf(d) === i),
      'the group labels, as the Deadline cells read');
  });

  test('field values: a Date in every stored shape reads as its date, a Boolean as true or false', async () => {
    const day = (v: unknown) => (fieldValue(v) as Date).toISOString();
    expect(day(formDate), formDate.date, 'the form\'s {_value, isUtc, date}');
    expect(day(oldDate), oldDate.date, 'the older {a, b, date}');
    expect(day({_value: Date.parse(formDate.date), isUtc: false}), formDate.date, 'only _value');
    expect(fieldValue({unknown: 1}) === undefined, true, 'an object of no known shape is no value');
    expect(fieldText(true), 'true', 'true');
    expect(fieldText(false), 'false', 'false');
    expect(fieldText('CCO'), 'CCO', 'text');
    expect(/^\d{4}\/\d{1,2}\/\d{1,2}$/.test(fieldText(formDate)), true,
      `a date field reads as the table's days (2026/3/27); it reads "${fieldText(formDate)}"`);
    expect(dayText(null), '', 'no date, no text (never 1970/1/1)');
  });
});
