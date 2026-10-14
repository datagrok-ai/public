/** "My programs" on Home: a tile per pinned program the user can see, across portfolios. A tile
 * says where the program stands on its track, how ready it is for its next gate, what needs a
 * look and whether new data landed. A click puts the program in the context panel; a double-click opens
 * its page in the Project Portal. The widget reads
 * again whenever Home comes back, so a pin made elsewhere shows without a reload. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {divV, span} from '@datagrok-libraries/u2';
import {pinnedTiles} from '../pinned-programs';
import {NOTHING_PINNED, ProgramTile} from '../model/pinning';
import {installStyles} from './styles';
import {openProgramPage, programCard} from './program-card';
import {showProgramInPanel} from './program-panel';

/** Kept for callers from before the card had a module of its own. */
export {openProgramPage};

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

/** A tile: the program's card. A click selects it and puts the program in the context panel, as the
 * platform does with what is clicked; a double-click opens its page. */
function tile(t: ProgramTile, select: (el: HTMLElement) => void): HTMLElement {
  const el = programCard(t, `pinned-tile-${slug(t.code)}`);
  el.title += '\nDouble-click to open its page';
  // Home takes a click as its own: stopped here, and the page opened after it is done with
  el.onclick = (e) => {
    e.stopPropagation();
    select(el);
    void showProgramInPanel(t.programId).catch((err) => grok.shell.error(err?.message ?? String(err)));
  };
  el.ondblclick = (e) => {
    e.stopPropagation();
    setTimeout(() => void openProgramPage(t).catch((err) => grok.shell.error(err?.message ?? String(err))), 50);
  };
  return el;
}

export class MyProgramsWidget extends DG.Widget {
  private readonly body = named(divV([], 'pp-tiles'), 'my-programs');
  private loading = false;

  constructor() {
    // two columns of the Home grid: room for tiles side by side
    super(divV([], 'pp-my-programs'));
    this.root.classList.add('widget-wide');
    installStyles();
    this.root.append(this.body);
    void this.refresh();
    // Home is kept while other views open; read again each time it comes back
    this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => setTimeout(() => {
      if (this.root.isConnected)
        void this.refresh();
    }, 50)));
  }

  async refresh(): Promise<void> {
    if (this.loading)
      return;
    this.loading = true;
    try {
      const tiles = await pinnedTiles();
      const select = (el: HTMLElement) => {
        for (const x of Array.from(this.body.querySelectorAll('.pp-selected')))
          x.classList.remove('pp-selected');
        el.classList.add('pp-selected');
      };
      this.body.replaceChildren(...(tiles.length ? tiles.map((t) => tile(t, select)) :
        [named(span(NOTHING_PINNED, 'pp-muted'), 'no-pinned-programs')]));
    } catch (e: any) {
      this.body.replaceChildren(named(span(`Pinned programs could not be read: ${e?.message ?? e}`, 'pp-reason'), 'pinned-error'));
    } finally {
      this.loading = false;
    }
  }
}
