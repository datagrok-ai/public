/** A program's platform record as the platform shows it where favorites are listed, the Workspace
 * on Home (Power Pack) first. The favorite comes back as a bare record (entity type
 * `EntityRecord`, named `portal.program:<portfolio id>|<code>`) that no platform meta claims, so
 * without this handler the Workspace can only show it in the context panel. With it, the record
 * gets the program's icon, a card with its portfolio in the context panel, and Open opens the
 * program's page in the Project Portal. Registered when the package loads (`portalInit`); the
 * "My programs" widget on Home loads the package. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {isProgramRecordName, programAddress, PROGRAM_RECORD} from './model/pinning';
import {portalAppView} from './views/portal-view';
import {openProgramPage} from './views/my-programs-widget';
import {portalIcon} from './views/icons';

/** The program and portfolio behind a record; null when the user cannot read them. */
async function programOfRecord(x: DG.Entity): Promise<{program: ProgramRow, portfolio: PortfolioRow} | null> {
  const program = await portalDb.programs.get(x.id).catch(() => null);
  const portfolio = program ? await portalDb.portfolios.get(program.portfolio_id).catch(() => null) : null;
  return program && portfolio ? {program, portfolio} : null;
}

export class ProgramRecordHandler extends DG.ObjectHandler<DG.Entity> {
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

  /** The context panel: the program, its portfolio, and a link that opens its page. */
  renderProperties(x: DG.Entity): HTMLElement {
    const host = ui.divV([ui.divText(this.getCaption(x))]);
    host.dataset.u2Name = 'program-record-card';
    void programOfRecord(x).then((r) => {
      if (!r) {
        host.append(ui.divText('Not readable in the Project Portal'));
        return;
      }
      const open = ui.link('Open in the Project Portal', () => void openProgramPage({portfolioCode: r.portfolio.code,
        programId: r.program.id, address: programAddress(r.portfolio.code, r.program.code)}));
      open.dataset.u2Name = 'program-record-open';
      host.replaceChildren(ui.divText(r.program.name || r.program.code), ui.divText(`${r.portfolio.name} · ${r.program.code}`), open);
    });
    return host;
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

/** Registers the handler once per session. */
export function registerProgramRecordHandler(): void {
  if (registered)
    return;
  registered = true;
  DG.ObjectHandler.register(new ProgramRecordHandler());
}
