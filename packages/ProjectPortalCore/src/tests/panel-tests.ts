import * as grok from 'datagrok-api/grok';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {MORE_SECTIONS, PAGE_PANELS, panelLayout, PersonaConfig, validatePersonas} from '../model/personas';
import {savedMode, saveMode, savedPersona, savePersona} from '../personas';
import {possessive, PortalView, startingPortfolio} from '../views/portal-view';
import {programPanel} from '../views/program-panel';
import {loadExampleDesignSets} from '../design-sets';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

const config: PersonaConfig = {panels: {lead: ['gate', 'approvals', 'feed'], biologist: ['dmta', 'compounds']}};

category('Portal: persona panel rules', () => {
  test('a composed persona: its panels in its order in both modes, the rest behind More', async () => {
    const l = panelLayout('lead', config, 'advanced', false);
    expect(`${l.composed} ${l.shown.join()}`, 'true gate,approvals,feed');
    expect(l.more.length, PAGE_PANELS.length - 3);
    expect(panelLayout('lead', config, 'simple', false).shown.join(), 'gate,approvals,feed', 'the mode does not hide them');
    const open = panelLayout('lead', config, 'simple', true);
    expect(open.shown.slice(0, 4).join(), 'gate,approvals,feed,quick_actions', 'More open: the rest after, in the page\'s order');
    expect(open.shown.length, PAGE_PANELS.length);
  });

  test('a persona the portfolio leaves out keeps the page; the mode decides as before', async () => {
    const a = panelLayout('chemist', config, 'advanced', false);
    expect(`${a.composed} ${a.shown.join() === PAGE_PANELS.join()}`, 'false true');
    const s = panelLayout('chemist', config, 'simple', false);
    expect(s.shown.some((x) => (MORE_SECTIONS as string[]).includes(x)), false);
    expect(s.more.join(), MORE_SECTIONS.join());
    expect(panelLayout(null, config, 'advanced', false).composed, false);
  });

  test('the config check: personas and panels it knows, once each', async () => {
    expect(validatePersonas(config).length, 0);
    expect(validatePersonas({panels: {lead: ['gate', 'gate'], boss: ['feed'], chemist: ['charts']}} as any).join('; '),
      'personas: lead lists a panel twice; personas: panels for "boss", which is not one of chemist, biologist, modeler, ' +
      'process_scientist, lead; personas: chemist panel "charts" is not one of ' + PAGE_PANELS.join(', '));
  });
});

category('Portal: persona panels', () => {
  let kin1: ProgramRow;
  let mode: ReturnType<typeof savedMode>;
  let persona: ReturnType<typeof savedPersona>;

  before(async () => {
    await loadAllPortfolios(_package);
    await loadExampleDesignSets(_package);
    const smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`, limit: ROW_LIMIT}))[0];
    mode = savedMode();
    persona = savedPersona();
  });

  after(async () => {
    saveMode(mode ?? 'advanced');
    savePersona(persona ?? null);
    await grok.userSettings.flush();
  });

  test('EX-SMD composes the biologist\'s page: DMTA, observations, compounds, activity; the gate behind More', async () => {
    const smd = await portfolioByCode('EX-SMD');
    expect(JSON.stringify((smd.personas as PersonaConfig).panels?.biologist),
      '["quick_actions","data","dmta","observations","compounds","feed"]', 'loaded from the file');
    saveMode('advanced');
    savePersona('biologist');
    const v = new PortalView();
    await v.load('EX-SMD');
    await v.showProgram('EX-SMD', kin1.id);
    const q = (name: string) => v.detailRoot.querySelector(`[data-u2-name="${name}"]`) as HTMLElement | null;
    expect(!!q('panels-composed') && !!q('dmta-heading') && !!q('observations-heading') && !!q('feed-heading'), true);
    expect(!q('gate-heading') && !q('saved-runs-heading') && !q('design-sets-heading'), true, 'not listed: behind More');
    const order = ['quick-actions-heading', 'dmta-heading', 'observations-heading', 'feed-heading']
      .map((n) => Array.from(v.detailRoot.querySelectorAll('[data-u2-name]')).indexOf(q(n)!));
    expect(order.every((x, i) => i === 0 || x > order[i - 1]), true, `in the listed order: ${order.join()}`);
    expect(q('more-button')?.textContent?.startsWith('More ('), true);
    q('more-button')!.click();
    expect(!!q('gate-heading') && q('more-button')?.textContent === 'Less', true);
    q('more-button')!.click();
    v.setPersona('chemist');
    expect(!q('panels-composed') && !!q('gate-heading') && !!q('saved-runs-heading') && !q('more-button'), true,
      'the chemist keeps the whole page');
  });
});

category('Portal: the board chooses, the program page works, the context panel acts', () => {
  test('a card puts the program in the context panel; opening it swaps the board for its page, and back', async () => {
    await loadAllPortfolios(_package);
    const smd = await portfolioByCode('EX-SMD');
    const kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`,
      limit: ROW_LIMIT}))[0];
    const v = new PortalView();
    await v.load('EX-SMD');
    const q = (name: string) => v.root.querySelector(`[data-u2-name="${name}"]`) as HTMLElement | null;
    const shown = (name: string) => q(name)!.style.display !== 'none';
    await v.previewProgram('EX-SMD', kin1.id);
    const o = grok.shell.o;
    expect(o?.id ?? o?.dataset?.u2Name, o?.id ? kin1.id : 'program-panel', 'the program\'s record (or its panel) is in the context panel');
    expect(shown('stage-board') && !shown('program-page'), true, 'the board stays');
    expect(q('program-exs-kin1')!.classList.contains('pp-selected'), true);
    await v.showProgram('EX-SMD', kin1.id);
    expect(!shown('stage-board') && shown('program-page'), true, 'the page takes the board\'s place');
    expect(v.path.value.endsWith('/EX-SMD/EXS-KIN1'), true);
    expect(!!q('gate-heading') && !!q('board-link'), true);
    expect(!q('program-actions'), true, 'the program\'s actions are in the context panel, not under its title');
    q('board-link')!.click();
    expect(shown('stage-board') && !shown('program-page'), true, 'back to the board');
    expect(v.path.value.endsWith('/EX-SMD'), true);
    expect(!q('gate-heading'), true, 'the page is taken down with its scope');
  });

  test('the program panel: its card, the favorite star, Open and the actions', async () => {
    await loadAllPortfolios(_package);
    const smd = await portfolioByCode('EX-SMD');
    const kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`,
      limit: ROW_LIMIT}))[0];
    const panel = programPanel(kin1.id);
    const q = (name: string) => panel.querySelector(`[data-u2-name="${name}"]`) as HTMLElement | null;
    for (let i = 0; i < 100 && !q('program-actions'); i++)
      await new Promise((r) => setTimeout(r, 100));
    expect(!!q('program-panel-card') && !!q('pin-button') && !!q('open-program-button'), true);
    expect(/to favorites|from favorites/.test(q('pin-button')!.textContent ?? ''), true);
    for (const name of ['follow-button', 'program-provenance-button', 'ask-sme-button', 'program-data-package-button'])
      expect(!!q(name), true, name);
  });
});

category('Portal: which portfolio the app opens on', () => {
  test('the one portfolio opens straight away; several leave the choice to the landing page', async () => {
    expect(startingPortfolio(['EX-SMD']), 'EX-SMD');
    expect(startingPortfolio(['EX-SMD', 'EX-BIOC']), null);
    expect(startingPortfolio([]), null);
  });

  test('a portfolio asked for by the address opens if the user can see it', async () => {
    expect(startingPortfolio(['EX-SMD', 'EX-BIOC'], 'EX-BIOC'), 'EX-BIOC');
    expect(startingPortfolio(['EX-SMD', 'EX-BIOC'], 'EX-BIOD'), null, 'one the user cannot see: the landing page');
    expect(startingPortfolio(['EX-SMD'], 'EX-BIOD'), 'EX-SMD', 'one the user cannot see, with one they can: that one');
  });

  test('a group\'s favorites: "Viewers\'", "Team\'s"', async () => {
    expect(possessive('Portal EX-BIOC Viewers'), 'Portal EX-BIOC Viewers\'');
    expect(possessive('Team'), 'Team\'s');
  });
});
