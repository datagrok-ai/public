/** A program in the platform's context panel: its card, its favorite star and what can be done with
 * it. The program's record shows this wherever the platform puts the record in the context panel: a
 * card clicked on the board, the program's page, a row of the portfolio table, a tile or a favorite on
 * Home. The star and the group pin are the platform's favorites (`DG.Favorites`), as everywhere. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {button, divH, divV, icon, span} from '@datagrok-libraries/u2';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {isApprover, programContext} from '../portal';
import {follow, myFollows, unfollow} from '../new-data';
import {pinnableRecord, pinStateOf, setPinned} from '../pinning';
import {buildPortfolioWorkspace, existingSpaces, summarizeWorkspace} from '../workspace';
import {programTile} from '../pinned-programs';
import {programAddress} from '../model/pinning';
import {smeContextFrom} from '../model/discussion';
import {openProvenanceDialog} from './provenance-dialog';
import {openDataPackageDialog} from './data-package-dialog';
import {openAskSmeDialog} from './discussion-dialog';
import {openPortal, possessive} from './portal-view';
import {openProgramPage, programCard, showOnBoard} from './program-card';
import {linkTo} from './u2-controls';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const failed = (e: any): void => {
  grok.shell.error(e?.message ?? String(e));
};

/** Puts the program in the context panel: its record, which the package's handler draws with this
 * panel (and the platform keeps as the current object); the panel itself where the platform has no
 * record for it. */
export async function showProgramInPanel(programId: string): Promise<void> {
  const record = await grok.dapi.entities.find(programId).catch(() => null);
  grok.shell.windows.showContextPanel = true;
  grok.shell.o = record ?? programPanel(programId);
}

/** The panel of [programId]: drawn at once, filled as its facts are read. */
export function programPanel(programId: string): HTMLElement {
  installStyles();
  const host = named(divV([span('Reading the program…', 'pp-muted')], 'pp-program-panel'), 'program-panel');
  void fill(host, programId).catch((e) => host.replaceChildren(named(span(`The program could not be read: ${e?.message ?? e}`,
    'pp-reason'), 'program-panel-error')));
  return host;
}

async function fill(host: HTMLElement, programId: string): Promise<void> {
  const program = await portalDb.programs.get(programId).catch(() => null);
  const portfolio = program ? await portalDb.portfolios.get(program.portfolio_id).catch(() => null) : null;
  if (!program || !portfolio) {
    host.replaceChildren(named(span('Not readable in the Project Portal', 'pp-muted'), 'program-panel-unreadable'));
    return;
  }
  const noun = (portfolio.program_noun ?? 'program').toLowerCase();
  const [tile, follows, held, pinReason, approver, live] = await Promise.all([
    programTile(programId).catch(() => null),
    myFollows(programId).catch(() => []),
    pinStateOf(program).catch(() => ({mine: false, groups: [] as string[]})),
    pinnableRecord(program).then((r) => r.reason).catch(() => null),
    isApprover(portfolio).catch(() => false),
    program.space_id ? existingSpaces([program.space_id]).then((s) => s.has(program.space_id!)).catch(() => false) :
      Promise.resolve(false),
  ]);
  const viewers = approver && portfolio.viewers ? await grok.dapi.groups.find(portfolio.viewers).catch(() => null) : null;
  // drawn again after a change, and the page too when it shows this program
  const changed = async () => {
    await openPortal()?.portal.refreshProgram(programId);
    if (host.isConnected)
      await fill(host, programId);
  };
  const address = {portfolioCode: portfolio.code, programId, address: programAddress(portfolio.code, program.code)};

  // the favorite star, as the platform's panels have it
  const mine = held.mine;
  const star = named(divH([icon('star', {variant: mine ? 'solid' : 'regular'}),
    span(mine ? 'Remove from favorites' : 'Add to favorites')], `pp-star${mine ? ' pp-star-on' : ''}`), 'pin-button');
  star.title = pinReason && !mine ? `Cannot add to favorites: ${pinReason}` :
    mine ? 'Take it out of your favorites and your Workspace, and its tile off Home' :
      'Put it in your favorites and your Workspace on Home, with a tile in My programs';
  if (!pinReason || mine) {
    star.onclick = () => void setPinned(portfolio, program, !mine, null)
      .then((line) => grok.shell.info(line)).then(changed).catch(failed);
  } else
    star.classList.add('pp-star-off');

  const shown = openPortal()?.portal;
  const onPage = !!shown?.showsPage(programId);
  const open = named(button(`Open ${noun}`, () => void openProgramPage(address).catch(failed), {primary: true}),
    'open-program-button');
  open.title = `This ${noun}'s page: gate, documents, runs, discussion and the rest`;

  const action = (text: string, name: string, title: string, run: () => unknown) =>
    linkTo(text, name, () => void Promise.resolve(run()).catch(failed), {title});
  const followed = follows.some((f) => !f.work_unit_code);
  const actions: HTMLElement[] = [
    ...(onPage || shown?.showsBoard(portfolio.code) ? [] : [action('Show in portfolio', 'show-in-portfolio-link',
      `Its place among the ${portfolio.name} programs, on the portfolio's board`, () => showOnBoard(address))]),
    action(followed ? 'Unfollow' : 'Follow', 'follow-button', followed ? 'Stop being notified of new data in this program' :
      'Be notified, in the platform\'s notifications, when new results land in this program\'s data', async () => {
      if (followed)
        await unfollow(program);
      else
        await follow(portfolio, program);
      grok.shell.info(`${program.name}: ${followed ? 'no longer followed' : 'followed; new results will show in your notifications'}`);
      await changed();
    }),
    ...(live ? [action('Open workspace', 'open-workspace-button', 'The program\'s folder, as the Browse tree opens it',
      () => openWorkspace(program))] : approver && portfolio.workspace_root ? [action('Build its folder', 'build-folder-button',
      'Create this program\'s folder and its workstream folders in the portfolio\'s workspace', async () => {
        grok.shell.info(summarizeWorkspace(await buildPortfolioWorkspace(portfolio.code, false, program.id)));
        await changed();
      })] : [named(span(program.space_id ? 'The workspace folder is gone' : 'No workspace yet', 'pp-muted'), 'no-workspace')]),
    ...(viewers ? [groupPin(portfolio, program, viewers, held.groups, pinReason, changed)] : []),
    action('Provenance', 'program-provenance-button', 'What the program\'s records rest on, and what uses them',
      () => openProvenanceDialog('program', program.id, {openProgram: async (id) => {
        const p = await portalDb.programs.get(id);
        await openProgramPage({portfolioCode: (await portalDb.portfolios.get(p.portfolio_id)).code, programId: id,
          address: programAddress((await portalDb.portfolios.get(p.portfolio_id)).code, p.code)});
      }})),
    action('Ask SME…', 'ask-sme-button', 'Ask an expert about this program: the portfolio, stage, next gate and readiness go with the question',
      async () => {
        const ctx = await programContext(program.id);
        openAskSmeDialog({portfolio, program, context: smeContextFrom({...ctx, anchor: null, anchorUrl: null,
          appUrl: `${window.location.origin}${address.address}`}), onAsked: changed});
      }),
    action('Data package', 'program-data-package-button',
      'What the program\'s gate case rests on, one line each, with status and source; exports to CSV',
      () => openDataPackageDialog(program.id)),
  ];
  host.replaceChildren(
    tile ? programCard(tile, 'program-panel-card', [star]) :
      divH([span(program.name, 'pp-card-title pp-grow'), star], 'pp-card-meta'),
    ...(onPage ? [] : [divH([open], 'pp-actions')]),
    span('Actions', 'pp-heading'),
    named(divV(actions, 'pp-panel-actions'), 'program-actions'),
  );
}

/** For an approver: the program in (or out of) the favorites of the portfolio's viewers group. */
function groupPin(portfolio: PortfolioRow, program: ProgramRow, g: DG.Group, groups: string[], reason: string | null,
  changed: () => Promise<void>): HTMLElement {
  const gName = g.friendlyName || g.name;
  const on = groups.includes(gName);
  return linkTo(on ? `Remove from ${possessive(gName)} favorites` : `Add to ${possessive(gName)} favorites`, 'pin-group-button',
    () => {
      if (reason && !on)
        return grok.shell.warning(`Cannot add to favorites: ${reason}`);
      void setPinned(portfolio, program, !on, g).then((line) => grok.shell.info(line)).then(changed).catch(failed);
    }, {title: `${on ? 'Take it out of' : 'Put it in'} the favorites and Workspace of everyone in ${gName}`});
}

/** Opens the program's folder the way the Browse tree does. */
async function openWorkspace(program: ProgramRow): Promise<void> {
  const space = await grok.dapi.spaces.find(program.space_id!).catch(() => null);
  if (!space)
    return grok.shell.warning('This program\'s workspace folder is gone; rebuild the portfolio\'s workspace');
  grok.shell.addView(DG.SpaceView.forProject(space));
}
