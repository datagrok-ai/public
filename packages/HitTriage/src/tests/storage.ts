/* eslint-disable max-len */
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {CampaignGrouping, HDCampaignsCollapsedLSKey, HDCampaignsGroupingLSKey, HDCampaignsTableSortingLSKey,
  HDCampaignTableColumnsLSKey, HDDefaultCampaignTableColumns, HTCampaignsCollapsedLSKey, HTCampaignsGroupingLSKey,
  HTCampaignsTableSortingLSKey, HTCampaignTableColumnsLSKey, HTDefaultCampaignTableColumns, HTFunctionOrderingLSKey} from '../app/consts';
import {getSavedFunctionOrdering, setSavedFunctionOrdering} from '../app/utils';
import {catalogPersistence} from '../app/shared/campaign-catalog';
import {expectSame} from './fixtures';

/** Runs [body] with [keys] cleared, then puts the viewer's own values back. */
async function withKeys(keys: string[], body: () => Promise<void> | void): Promise<void> {
  const saved = keys.map((k) => localStorage.getItem(k));
  try {
    for (const k of keys)
      localStorage.removeItem(k);
    await body();
  } finally {
    for (const [i, k] of keys.entries())
      saved[i] === null ? localStorage.removeItem(k) : localStorage.setItem(k, saved[i]!);
  }
}

const HD_KEYS = [HDCampaignsGroupingLSKey, HDCampaignTableColumnsLSKey, HDCampaignsTableSortingLSKey, HDCampaignsCollapsedLSKey];
const HT_KEYS = [HTCampaignsGroupingLSKey, HTCampaignTableColumnsLSKey, HTCampaignsTableSortingLSKey, HTCampaignsCollapsedLSKey];
const state = (patch: object) => ({columns: [], sort: '', groupBy: null, collapsed: {}, ...patch});
const groupKey = (key: string) => `Accordion:Hit-Design-campaigns-group-${key}`;

category('Storage', () => {
  test('the grouping keys hold the grouping as plain text', () => withKeys([...HD_KEYS, ...HT_KEYS], () => {
    const hd = catalogPersistence('HD');
    const ht = catalogPersistence('HT');
    expect(hd.load().groupBy, null, 'Hit Design default');
    expect(ht.load().groupBy, null, 'Hit Triage default');
    hd.save(state({groupBy: CampaignGrouping.LastModifiedUser}));
    expect(localStorage.getItem('HDCampaignsGrouping'), 'Last Modified User', 'HDCampaignsGrouping');
    ht.save(state({groupBy: 'campaignFields.Chemist'}));
    expect(localStorage.getItem('HTCampaignsGrouping'), 'campaignFields.Chemist', 'HTCampaignsGrouping');
    localStorage.setItem('HDCampaignsGrouping', 'campaignFields.Target MW');
    expect(hd.load().groupBy, 'campaignFields.Target MW', 'a stored field grouping reads back');
    expect(ht.load().groupBy, 'campaignFields.Chemist', 'the two apps keep separate keys');
    hd.save(state({}));
    expect(localStorage.getItem('HDCampaignsGrouping'), CampaignGrouping.None, 'no grouping is "None"');
    expect(hd.load().groupBy, null, '"None" reads back as no grouping');
  }));

  test('the column keys hold a JSON list of column names in display order', () => withKeys([...HD_KEYS, ...HT_KEYS], () => {
    const hd = catalogPersistence('HD');
    expectSame(HDDefaultCampaignTableColumns, ['Code', 'Created', 'Author', 'Last Modified by', 'Molecules', 'Status'], 'Hit Design default');
    expectSame(HTDefaultCampaignTableColumns, ['Created', 'Author', 'Total', 'Selected', 'Status'], 'Hit Triage default');
    expect(hd.load().columns === undefined, true, 'nothing kept: the defaults');
    hd.save(state({columns: ['Status', 'Code', 'campaignFields.Chemist']}));
    expect(localStorage.getItem('HDCampaignTableColumns'), '["Status","Code","campaignFields.Chemist"]', 'HDCampaignTableColumns');
    catalogPersistence('HT').save(state({columns: ['Selected', 'Total']}));
    expect(localStorage.getItem('HTCampaignTableColumns'), '["Selected","Total"]', 'HTCampaignTableColumns');
    expectSame(hd.load().columns, ['Status', 'Code', 'campaignFields.Chemist'], 'read back in order');
    for (const bad of ['{not json', '[1,2]', '"Code"', '{"a":1}']) {
      localStorage.setItem('HDCampaignTableColumns', bad);
      expect(hd.load().columns === undefined, true, `${bad} falls back to the default`);
    }
  }));

  test('the sorting keys hold columnName and ascending, or nothing', () => withKeys([...HD_KEYS, ...HT_KEYS], () => {
    const hd = catalogPersistence('HD');
    const ht = catalogPersistence('HT');
    expect(hd.load().sort, '', 'Hit Design default');
    expect(ht.load().sort, '', 'Hit Triage default');
    hd.save(state({sort: '!Created'}));
    expect(localStorage.getItem('HDCampaignsTableSorting'), '{"columnName":"Created","ascending":false}', 'HDCampaignsTableSorting');
    localStorage.setItem('HTCampaignsTableSorting', '{"columnName":"campaignFields.Chemist","ascending":true}');
    expect(ht.load().sort, 'campaignFields.Chemist', 'HTCampaignsTableSorting read back');
    hd.save(state({sort: ''}));
    expect(localStorage.getItem('HDCampaignsTableSorting'), null, 'no sorting removes the key');
  }));

  test('a save writes only the choice that changed: the default columns are not frozen', () => withKeys(HD_KEYS, () => {
    const hd = catalogPersistence('HD');
    hd.load();
    hd.save(state({columns: HDDefaultCampaignTableColumns, sort: '!Created'}));
    expect(localStorage.getItem('HDCampaignsTableSorting'), '{"columnName":"Created","ascending":false}', 'the sort chosen is kept');
    expect(localStorage.getItem('HDCampaignTableColumns'), null, 'the default columns are not written');
    expect(localStorage.getItem('HDCampaignsGrouping'), null, 'nor the grouping nobody chose');
    expect(localStorage.getItem(HDCampaignsCollapsedLSKey), null, 'nor the folds');
  }));

  test('the compute tab order key holds order and hidden lists', () => withKeys([HTFunctionOrderingLSKey], () => {
    expectSame(getSavedFunctionOrdering(), {order: [], hidden: []}, 'default');
    setSavedFunctionOrdering({order: ['Chem:addChemPropertiesColumns', 'descriptors'], hidden: ['Chembl:chemblMolregno']});
    expectSame(JSON.parse(localStorage.getItem('HTFunctionOrderingLS')!),
      {order: ['Chem:addChemPropertiesColumns', 'descriptors'], hidden: ['Chembl:chemblMolregno']}, 'HTFunctionOrderingLS');
    localStorage.setItem('HTFunctionOrderingLS', '{"order":["a"]}');
    expectSame(getSavedFunctionOrdering(), {order: ['a'], hidden: []}, 'a partial value reads with an empty hidden list');
  }));

  test('folded groups are kept per grouping; the accordion\'s own folds are read once', () =>
    withKeys([...HD_KEYS, ...Object.keys(localStorage).filter((k) => k.startsWith(groupKey(''))), groupKey('ZZB3 A'),
      groupKey('ZZB3 B')], () => {
      const hd = catalogPersistence('HD');
      localStorage.setItem('HDCampaignsGrouping', CampaignGrouping.Status);
      localStorage.setItem(groupKey('ZZB3 A'), '{"ZZB3 A":false}');
      localStorage.setItem(groupKey('ZZB3 B'), '{"ZZB3 B":true}');
      expectSame(hd.load().collapsed, {Status: ['ZZB3 A']}, 'only a group the accordion kept closed starts folded');
      hd.save(state({groupBy: CampaignGrouping.Status, collapsed: {Status: ['ZZB3 B'], Template: ['KHL']}}));
      expect(localStorage.getItem(HDCampaignsCollapsedLSKey), '{"Status":["ZZB3 B"],"Template":["KHL"]}', 'HDCampaignsCollapsed');
      expectSame(hd.load().collapsed, {Status: ['ZZB3 B'], Template: ['KHL']}, 'kept folds win over the accordion\'s');
    }));
});
