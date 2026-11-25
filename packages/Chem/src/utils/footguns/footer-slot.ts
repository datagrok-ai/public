/**
 * The middle slot of a dialog's footer, where the foot-gun flags of what OK would write show before OK (the footer has
 * room in its middle, and the buttons must not jump around). The slot is made with the dialog and stays in its footer,
 * empty or not, one line high, as wide as the space the footer has between its left edge and its buttons and never as
 * wide as its text: so OK and Cancel stay where they are whatever it says. A long text is cut with an ellipsis; the
 * slot's tooltip (and its `aria-label`) holds all of it. A block disables OK (the dialog's own check: a disabled OK is
 * neither clicked nor run by Enter); a warn or a note leaves OK enabled. What the flags are and how they are worded is
 * checks.ts's.
 */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Flag, LEVEL_ORDER, Level, flagsText, worstLevel} from './checks';

/** The slot's height (CSS px): one line of the footer's text; the footer itself is 44 px. */
const SLOT_HEIGHT = 20;
/** The space between the slot and the buttons, and the dialog's left edge (CSS px). */
const SLOT_MARGIN = 8;
/** The slot's text: its size (CSS px) and its colour per level (Datagrok's palette variables). */
const SLOT_FONT_SIZE = 12;
const LEVEL_COLOUR: {[level in Level]: string} = {block: 'var(--red-3)', warn: 'var(--orange-2)', note: 'var(--grey-5)'};
/** The icon before the text, per level (Font Awesome names). */
const LEVEL_ICON: {[level in Level]: string} = {block: 'ban', warn: 'exclamation-triangle', note: 'info-circle'};
/** The slot's class: what the tests and the BDD bindings find it by. */
export const FOOTER_SLOT_CLASS = 'chem-footgun-slot';
/** The class of the words beside the Reaction Enumerator's Enumerate button (the same line, in its ribbon). */
export const ENUMERATE_WORDS_CLASS = 'chem-footgun-enumerate';
/** Those words' widest (CSS px): the ribbon sizes its panels by their content. */
const WORDS_MAX_WIDTH = 420;

/** The middle slot of one dialog's footer; made with no dialog, the same line standing alone (the words beside the
 * enumerator's Enumerate button), for its owner to place. */
export class FooterSlot {
  /** The slot element, in the footer from the start. */
  readonly root: HTMLDivElement;
  private readonly icon = ui.div([], {style: {flex: '0 0 auto', marginRight: '4px'}});
  private readonly label = ui.div([], {style: {flex: '1 1 0', minWidth: '0', overflow: 'hidden',
    textOverflow: 'ellipsis', whiteSpace: 'nowrap'}});
  private flags: Flag[] = [];
  // The dialog's OK, found as the slot is made (a dialog may rename it later: New Reaction's SAVE).
  private readonly ok: HTMLButtonElement | null;

  constructor(dialog: DG.Dialog | null, className = FOOTER_SLOT_CLASS) {
    this.root = ui.div([this.icon, this.label], className);
    // its size the footer's, never its text's: it grows into the free space from a zero basis, and its least width is
    // zero, so neither a text nor its absence changes where the buttons are or how wide the dialog is
    Object.assign(this.root.style, {display: 'flex', alignItems: 'center', flex: '1 1 0', minWidth: '0',
      height: `${SLOT_HEIGHT}px`, lineHeight: `${SLOT_HEIGHT}px`, margin: `0 ${SLOT_MARGIN}px`, overflow: 'hidden',
      fontSize: `${SLOT_FONT_SIZE}px`});
    if (dialog === null) {
      // standing alone (a ribbon): as wide as its text, up to a limit (its text from its own width, not from zero)
      Object.assign(this.root.style, {flex: '0 1 auto', maxWidth: `${WORDS_MAX_WIDTH}px`});
      this.label.style.flex = '0 1 auto';
      this.ok = null;
    } else {
      // a dialog sized to its content (most are: absolutely placed, shrink-to-fit) would grow to fit a long text and
      // move its buttons (as a package test found): the slot's own width is zero to its container's sizing, and it
      // takes only the space the footer has to spare
      Object.assign(this.root.style, {width: '0', contain: 'inline-size'});
      // the footer's buttons are laid out right to left (OK rightmost): the slot goes after the last of them in the
      // markup, so it shows left of Cancel, before any icon the footer keeps at its left edge
      this.ok = dialog.getButton('OK') ?? null;
      const bar = this.ok?.parentElement ?? dialog.root.querySelector('.d4-command-bar');
      const buttons = bar ? [...bar.querySelectorAll(':scope > button')] : [];
      if (buttons.length)
        buttons[buttons.length - 1].after(this.root);
      else
        bar?.append(this.root);
    }
    ui.tooltip.bind(this.root, () => this.flags.length ? flagsText(this.flags, '\n') : null);
    this.show([]);
  }

  /** The flags shown now. */
  get shown(): readonly Flag[] {
    return this.flags;
  }

  /** Whether a block is shown: OK is disabled. */
  get blocked(): boolean {
    return this.flags.some((f) => f.level === 'block');
  }

  /** Shows `flags` (the worst first; empty: the slot is blank), and enables OK unless one is a block or `okOff` (the
   * dialog says itself, once, why it has nothing to write: the reaction editor's refusal, RXN-06). */
  show(flags: Flag[], okOff = false): void {
    this.flags = [...flags].sort((a, b) => LEVEL_ORDER[a.level] - LEVEL_ORDER[b.level]);
    const level = worstLevel(this.flags);
    ui.empty(this.icon);
    if (level !== null)
      this.icon.append(ui.iconFA(LEVEL_ICON[level]));
    const text = flagsText(this.flags, ' · ');
    this.label.textContent = text;
    this.root.style.color = level !== null ? LEVEL_COLOUR[level] : '';
    this.root.dataset.level = level ?? '';
    this.root.setAttribute('aria-label', flagsText(this.flags, '\n'));
    if (this.ok)
      this.ok.disabled = this.blocked || okOff;
  }
}

/** The footer slot of `dialog`, made now (once per dialog: call it as the dialog is built, before OK is renamed). */
export function footerSlot(dialog: DG.Dialog): FooterSlot {
  return new FooterSlot(dialog);
}
