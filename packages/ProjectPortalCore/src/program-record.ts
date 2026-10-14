/** A program's platform record as the platform shows it wherever it lists or selects it: the
 * Workspace and Favorites on Home (Power Pack), Browse, and the Project Portal itself, which puts the
 * record in the context panel for the program clicked on the board or opened on its page. The
 * favorite comes back as a bare record (entity type `EntityRecord`, named
 * `portal.program:<portfolio id>|<code>`) that no platform meta claims; this handler gives it the
 * program's icon, its card (`renderCard`, the Workspace's preview), its panel in the context panel
 * (card, favorite star, actions), and Open, which opens the program's page. Its context menu gets
 * Open and Show in portfolio. Registered when the package loads (`portalInit`). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {isProgramRecordName, programAddress, PROGRAM_RECORD} from './model/pinning';
import {openPortal, portalAppView} from './views/portal-view';
import {pinnedByMeNow, pinStateOf, refreshMyPins, setPinned} from './pinning';
import {openProgramPage, programCard, showOnBoard} from './views/program-card';
import {programPanel} from './views/program-panel';
import {programTile} from './pinned-programs';
import {portalIcon} from './views/icons';
import {installStyles} from './views/styles';

/** The program and portfolio behind a record; null when the user cannot read them. */
async function programOfRecord(x: DG.Entity): Promise<{program: ProgramRow, portfolio: PortfolioRow} | null> {
  const program = await portalDb.programs.get(x.id).catch(() => null);
  const portfolio = program ? await portalDb.portfolios.get(program.portfolio_id).catch(() => null) : null;
  return program && portfolio ? {program, portfolio} : null;
}

export class ProgramRecordHandler extends DG.ObjectHandler<DG.Entity> {
  /** Power Pack leaves Group favorites (every group on the server) off the record's menu: the user's own
   * favorite is on the menu, and an approver pins a program for the portfolio's viewers from its panel. */
  readonly ownsFavorites = true;

  get type(): string {
    return `${PROGRAM_RECORD} record`;
  }

  isApplicable(x: any): boolean {
    return x instanceof DG.Entity && x.entityType === 'EntityRecord' && isProgramRecordName(x.name);
  }

  getCaption(x: DG.Entity): string {
    return x.friendlyName || x.name;
  }

  renderIcon(x: DG.Entity): HTMLElement {
    return portalIcon('program', `Program ${this.getCaption(x)}`);
  }

  /** The icon and the name: lists that show only markup (Favorites on Home) get the icon too. */
  renderMarkup(x: DG.Entity): HTMLElement {
    return ui.divH([portalIcon('program'), ui.span([this.getCaption(x)])], {style: {alignItems: 'center'}});
  }

  renderTooltip(x: DG.Entity): HTMLElement {
    return ui.divText(`Program ${this.getCaption(x)}: double-click or Open to open its page in the Project Portal`);
  }

  /** The program's card: the Workspace on Home previews a selected favorite with it. */
  renderCard(x: DG.Entity): HTMLElement {
    installStyles();
    const host = ui.divV([ui.divText(this.getCaption(x))]);
    host.dataset.u2Name = 'program-record-card';
    void programTile(x.id).then((t) => host.replaceChildren(t ? programCard(t, 'program-record-tile') :
      ui.divText('Not readable in the Project Portal')));
    return host;
  }

  /** The context panel: the program's card, its favorite star and its actions. */
  renderProperties(x: DG.Entity): HTMLElement {
    return programPanel(x.id);
  }

  /** Open: the program's page in a Project Portal view. */
  async renderPreview(x: DG.Entity): Promise<DG.View> {
    const r = await programOfRecord(x);
    if (!r)
      throw new Error(`${this.getCaption(x)} is not readable in the Project Portal`);
    return await portalAppView(r.portfolio.code, r.program.code) as DG.View;
  }
}

let registered = false;

/** The program's page and its board place, by its record. */
async function addressOf(x: DG.Entity): Promise<{portfolioCode: string, programId: string, address: string}> {
  const r = await programOfRecord(x);
  if (!r)
    throw new Error(`${x.friendlyName || x.name} is not readable in the Project Portal`);
  return {portfolioCode: r.portfolio.code, programId: r.program.id, address: programAddress(r.portfolio.code, r.program.code)};
}

/** Registers the handler once per session, and the program record's context menu wherever the platform
 * shows it (the Workspace and Favorites on Home, Browse): Open and Show in portfolio first, then Add to /
 * Remove from favorites in place of the platform's Add To Favorites (which does not know the record is a
 * favorite already), then the platform's Create URL alias. Power Pack's Group favorites is left off
 * (`ownsFavorites`). */
export function registerProgramRecordHandler(): void {
  if (registered)
    return;
  registered = true;
  const handler = new ProgramRecordHandler();
  DG.ObjectHandler.register(handler);
  const failed = (e: any) => grok.shell.error(e?.message ?? String(e));
  grok.events.onContextMenu.subscribe((e: any) => {
    const x = e?.args?.item;
    if (!handler.isApplicable(x))
      return;
    const menu: DG.Menu = e.args.menu;
    // taken out, to come back below ours
    const take = (text: string) => {
      const item = menu.find(text);
      if (item)
        menu.remove(text);
      return item;
    };
    const alias = take('Create URL alias...');
    take('Add To Favorites');
    take('Remove From Favorites');
    menu.item('Open', () => void addressOf(x).then(openProgramPage).catch(failed));
    menu.item('Show in portfolio', () => void addressOf(x).then(showOnBoard).catch(failed));
    // favorites not read yet (the menu right after the platform starts): its state is read on the click
    let known = pinnedByMeNow([x.id]);
    const toggle = async () => {
      const r = await programOfRecord(x);
      if (!r)
        throw new Error(`${x.friendlyName || x.name} is not readable in the Project Portal`);
      const pinned = known ?? (await pinStateOf(r.program)).mine;
      grok.shell.info(await setPinned(r.portfolio, r.program, !pinned, null));
      await openPortal()?.portal.refreshProgram(r.program.id);
      // the program's panel, if it is the one on show, is drawn again with its star
      if ((grok.shell.o as any)?.id === x.id) {
        grok.shell.o = null;
        grok.shell.o = x;
      }
    };
    const UNKNOWN = 'Add to or remove from favorites';
    menu.item(known ? 'Remove from favorites' : known === false ? 'Add to favorites' : UNKNOWN,
      () => void toggle().catch(failed));
    const ALIAS = 'Create URL alias...';
    if (alias)
      menu.item(ALIAS, () => alias.click());
    // favorites not read yet (the menu right after the platform starts): the item says which once read
    if (known === null) {
      void refreshMyPins().then(() => {
        known = pinnedByMeNow([x.id]);
        if (known === null)
          return;
        menu.remove(UNKNOWN);
        if (alias)
          menu.remove(ALIAS);
        menu.item(known ? 'Remove from favorites' : 'Add to favorites', () => void toggle().catch(failed));
        if (alias)
          menu.item(ALIAS, () => alias.click());
      });
      return;
    }
    // for the next menu: a favorite changed elsewhere (the Workspace's own x) shows there
    void refreshMyPins();
  });
  // the package can load before the session is ready: read until it answers, a few times
  const firstRead = (left: number): void => void refreshMyPins().then(() => {
    if (pinnedByMeNow([]) === null && left > 0)
      setTimeout(() => firstRead(left - 1), 3000);
  });
  firstRead(5);
}
