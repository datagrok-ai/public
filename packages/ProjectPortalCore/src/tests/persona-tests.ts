import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {
  asMode, asPersona, DEFAULT_CARDS, defaultPersona, MORE_SECTIONS, PersonaConfig, quickCards, sectionShown, validatePersonas,
} from '../model/personas';
import {personaIn, PREFERENCES, runFunctionCard, savedMode, saveMode, savedPersona, savePersona} from '../personas';
import {PortalView} from '../views/portal-view';
import {_package} from '../package-test';

const config: PersonaConfig = {
  default: 'biologist',
  groups: [{group: 'Leads', persona: 'lead'}, {group: 'Chemists', persona: 'chemist'}],
  cards: {modeler: ['save_run', {name: 'Portfolio as a table', function: 'ProjectPortalCore:portfolioTable'}, 'save_run']},
};

category('Portal: persona rules', () => {
  test('the persona of a user who has not chosen: first mapped group, else the default, else by setting', async () => {
    expect(defaultPersona(config, ['All users', 'Chemists', 'Leads'], 'discovery'), 'lead', 'the first mapping wins');
    expect(defaultPersona(config, ['chemists'], 'discovery'), 'chemist', 'group names ignore case');
    expect(defaultPersona(config, ['All users'], 'development'), 'biologist');
    expect(defaultPersona(undefined, [], 'discovery'), 'chemist');
    expect(defaultPersona(null, [], 'development'), 'process_scientist');
    expect(`${asPersona('modeler')} ${asPersona('pharmacist')} ${asMode('simple')} ${asMode('')}`, 'modeler null simple null');
  });

  test('cards: the defaults, the portfolio\'s own, function cards, once each, and why one is off', async () => {
    expect(quickCards('chemist', undefined).map((c) => c.key).join(), DEFAULT_CARDS.chemist.join());
    const modeler = quickCards('modeler', config);
    expect(modeler.map((c) => c.key).join(), 'save_run,fn-ProjectPortalCore:portfolioTable');
    expect(`${modeler[1].label} ${modeler[1].function} ${modeler[1].action}`, 'Portfolio as a table ProjectPortalCore:portfolioTable null');
    expect(quickCards('lead', config).map((c) => c.key).join(), DEFAULT_CARDS.lead.join(), 'a persona left out keeps the defaults');
    const off = quickCards('modeler', undefined, {choose_run: 'No criterion of this gate rests on a run'},
      {open_data: {label: 'Open SAR table'}});
    expect(off.find((c) => c.key === 'choose_run')!.unavailable, 'No criterion of this gate rests on a run');
    expect(off.find((c) => c.key === 'open_data')!.label, 'Open SAR table');
    expect(off.find((c) => c.key === 'save_run')!.unavailable, null);
  });

  test('a persona config is checked: personas, action codes, function cards', async () => {
    expect(validatePersonas(config).join('; '), '');
    expect(validatePersonas(undefined).length, 0);
    const errors = validatePersonas({default: 'pharmacist' as any, groups: [{group: ' ', persona: 'lead'}],
      cards: {chemist: ['open_sar' as any, {name: 'X', function: 'nopackage'}], nurse: []} as any});
    expect(errors.length, 5, errors.join('; '));
    expect(errors.some((e) => e.includes('card "open_sar" is not one of')), true);
    expect(errors.some((e) => e.includes('needs a name and a function as Package:Function')), true);
  });

  test('Simple keeps analytics, runs, snapshots and history behind More; Advanced shows everything', async () => {
    expect(MORE_SECTIONS.join(), 'recommendation,attention,saved_runs,roles,snapshots,feed,history');
    expect(MORE_SECTIONS.every((s) => sectionShown('advanced', s, false)), true);
    expect(MORE_SECTIONS.some((s) => sectionShown('simple', s, false)), false);
    expect(MORE_SECTIONS.every((s) => sectionShown('simple', s, true)), true);
  });

  test('the example files map their approvers to the lead persona', async () => {
    const files = await readPortfolioFiles(_package);
    for (const code of ['EX-SMD', 'EX-BIOC']) {
      const c = files.find((f) => f.code === code)!;
      expect(validatePortfolio(c).join('; '), '');
      expect(c.personas?.groups?.[0]?.persona, 'lead', code);
    }
    const smd = files.find((f) => f.code === 'EX-SMD')!;
    expect(validatePortfolio({...smd, personas: {default: 'nurse' as any}}).some((e) => e.includes('personas: default')), true);
  });
});

category('Portal: modes and quick actions', () => {
  let kin1: ProgramRow;
  let was: {[key: string]: string} | undefined;
  const opened: DG.ViewBase[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    const smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    was = grok.userSettings.get(PREFERENCES);
  });

  after(async () => {
    savePersona((was?.persona as any) ?? null);
    saveMode((was?.mode as any) ?? 'advanced');
    if (!was?.mode)
      grok.userSettings.delete(PREFERENCES, 'mode');   // there: just written
    await grok.userSettings.flush();
    for (const v of opened.splice(0))
      v.close();
  });

  test('mode and persona are remembered in the user\'s settings', async () => {
    saveMode('simple');
    savePersona('modeler');
    await grok.userSettings.flush();
    expect(`${savedMode()} ${savedPersona()} ${grok.userSettings.getValue(PREFERENCES, 'persona')}`, 'simple modeler modeler');
    savePersona(null);
    expect(savedPersona(), null, 'forgotten: the portfolio decides again');
    const smd = await portfolioByCode('EX-SMD');
    expect(JSON.stringify(smd.personas?.groups), '[{"group":"Portal EX-SMD Approvers","persona":"lead"}]', 'loaded');
    const p = await personaIn(smd);
    expect(p.chosen, false);
    expect(['lead', 'chemist'].includes(p.persona), true, p.persona);
  });

  test('the program page: persona cards, then Simple hides the heavier sections behind More', async () => {
    saveMode('advanced');
    savePersona('chemist');
    const v = new PortalView();
    await v.load('EX-SMD');
    await v.showProgram('EX-SMD', kin1.id);
    const q = (name: string) => v.root.querySelector(`[data-u2-name="${name}"]`) as HTMLElement | null;
    expect(q('quick-actions-heading')?.textContent, 'Quick actions · Chemist');
    expect(q('quick-open-data')?.textContent?.startsWith('Open SAR table'), true, q('quick-open-data')?.textContent ?? 'none');
    expect(!!q('quick-needs-attention') && !q('quick-save-run'), true);
    expect(!!q('saved-runs-heading') && !q('more-button'), true, 'Advanced shows everything');

    v.setPersona('modeler');
    expect(q('quick-actions-heading')?.textContent, 'Quick actions · Modeler');
    expect(!!q('quick-save-run') && !!q('quick-suggested-models') && !q('quick-needs-attention'), true);
    expect(q('quick-choose-run')?.title, 'No criterion of this gate rests on a run');
    expect(savedPersona(), 'modeler');

    v.setMode('simple');
    expect(!q('saved-runs-heading') && !q('transition-history') && !q('suggested-models-heading'), true);
    expect(q('more-button')?.textContent, `More (${MORE_SECTIONS.length})`);
    expect(!!q('gate-heading') && !!q('quick-actions'), true, 'the gate and the cards stay');
    q('more-button')!.click();
    expect(!!q('saved-runs-heading') && q('more-button')?.textContent === 'Less', true);
    expect(savedMode(), 'simple');
    v.setMode('advanced');
    expect(!!q('saved-runs-heading') && !q('more-button'), true);
  });

  test('a function card runs its function for the program and opens the table it returns', async () => {
    const before = Array.from(grok.shell.views).length;
    await runFunctionCard('ProjectPortalCore:portfolioTable', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    const views = Array.from(grok.shell.views);
    expect(views.length, before + 1);
    const tv = views[views.length - 1] as DG.TableView;
    opened.push(tv);
    expect(tv.dataFrame.rowCount > 0 && tv.dataFrame.col('program') !== null, true);
    let refused = '';
    await runFunctionCard('ProjectPortalCore:noSuchFunction', {portfolio: 'EX-SMD', program: 'EXS-KIN1'})
      .catch((e) => refused = e.message);
    expect(refused, 'ProjectPortalCore:noSuchFunction is not installed here');
  });
});
