/** Collapsible KNIME-style node groups — editor-level like annotations: the graph never
 *  changes, members stay real nodes with real connections. */

import {setTid} from '../utils/test-ids';

/** Title-bar height of the expanded frame (canvas px); kept in sync with `.ff-group-titlebar` CSS. */
export const GROUP_TITLE_H = 26;
/** Padding between the member bounding box and the frame border. */
export const GROUP_PAD = 14;
/** First boundary-dot center and step on the minimized card — the editor's wire-anchor
 *  math and the card's dot rendering both derive from these. */
export const GROUP_DOT_TOP = 15;
export const GROUP_DOT_STEP = 14;

export interface GroupDoc {
  id: string;
  title: string;
  description: string;
  memberIds: string[];
  minimized: boolean;
  /** Card anchor (canvas coords); while expanded tracks the frame's top-left. */
  pos: {x: number; y: number};
}

export class FlowGroup {
  id: string;
  title: string;
  description: string;
  memberIds: Set<string>;
  minimized: boolean;
  pos: {x: number; y: number};
  /** Frame size while expanded (derived from the member bbox by the editor). */
  frameSize: {w: number; h: number} = {w: 200, h: 120};

  readonly element: HTMLElement;
  readonly titleBar: HTMLElement;
  readonly caretEl: HTMLElement;
  readonly titleEl: HTMLElement;
  readonly statusEl: HTMLElement;
  readonly descEl: HTMLElement;
  private readonly inDots: HTMLElement;
  private readonly outDots: HTMLElement;

  constructor(opts: Partial<GroupDoc> = {}) {
    this.id = opts.id ?? `grp-${Math.random().toString(36).slice(2, 10)}`;
    this.title = opts.title ?? 'Group';
    this.description = opts.description ?? '';
    this.memberIds = new Set(opts.memberIds ?? []);
    this.minimized = opts.minimized === true;
    this.pos = opts.pos ?? {x: 0, y: 0};

    this.element = document.createElement('div');
    this.element.className = 'ff-group';
    this.element.dataset.groupId = this.id;
    setTid(this.element, 'group');

    this.titleBar = document.createElement('div');
    this.titleBar.className = 'ff-group-titlebar';

    this.statusEl = document.createElement('span');
    this.statusEl.className = 'ff-group-status';
    setTid(this.statusEl, 'group-status');
    this.titleBar.appendChild(this.statusEl);

    this.titleEl = document.createElement('div');
    this.titleEl.className = 'ff-group-title';
    setTid(this.titleEl, 'group-title');
    this.titleEl.textContent = this.title;
    // The title bar is the drag handle — renaming is a double-click, not an always-armed editable.
    this.titleEl.contentEditable = 'false';
    this.titleEl.spellcheck = false;
    this.titleEl.addEventListener('input', () => {
      this.title = this.titleEl.textContent ?? '';
    });
    this.titleEl.addEventListener('dblclick', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      this.startTitleEdit();
    });
    this.titleEl.addEventListener('blur', () => {
      this.titleEl.contentEditable = 'false';
    });
    this.titleBar.appendChild(this.titleEl);

    this.caretEl = document.createElement('span');
    this.caretEl.className = 'ff-group-caret';
    setTid(this.caretEl, 'group-caret');
    this.titleBar.appendChild(this.caretEl);

    this.element.appendChild(this.titleBar);

    this.descEl = document.createElement('div');
    this.descEl.className = 'ff-group-desc';
    setTid(this.descEl, 'group-desc');
    this.descEl.textContent = this.description;
    this.descEl.contentEditable = 'false';
    this.descEl.spellcheck = false;
    this.descEl.addEventListener('input', () => {
      this.description = this.descEl.textContent ?? '';
    });
    this.descEl.addEventListener('dblclick', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      this.startDescEdit();
    });
    this.element.appendChild(this.descEl);

    this.inDots = document.createElement('div');
    this.inDots.className = 'ff-group-dots ff-group-dots-in';
    this.element.appendChild(this.inDots);
    this.outDots = document.createElement('div');
    this.outDots.className = 'ff-group-dots ff-group-dots-out';
    this.element.appendChild(this.outDots);

    this.applyMode();
  }

  applyMode(): void {
    this.element.dataset.minimized = this.minimized ? 'true' : 'false';
    this.caretEl.textContent = this.minimized ? '▸' : '▾';
    this.caretEl.title = this.minimized ? 'Maximize' : 'Minimize';
    if (this.minimized) {
      this.element.style.left = `${this.pos.x}px`;
      this.element.style.top = `${this.pos.y}px`;
      this.element.style.width = '';
      this.element.style.height = '';
    } else {
      this.applyFrame();
    }
    // The card always shows the description; the frame hides an empty one until edited.
    this.descEl.style.display =
      this.description === '' && !this.descEl.dataset.editing ? 'none' : '';
  }

  applyFrame(): void {
    this.element.style.left = `${this.pos.x}px`;
    this.element.style.top = `${this.pos.y}px`;
    this.element.style.width = `${this.frameSize.w}px`;
    this.element.style.height = `${this.frameSize.h}px`;
  }

  applyCardPos(): void {
    this.element.style.left = `${this.pos.x}px`;
    this.element.style.top = `${this.pos.y}px`;
  }

  /** Aggregate member run status → title-bar dot (CSS keys off data-status). */
  setStatus(status: string): void {
    this.element.dataset.status = status;
    this.statusEl.dataset.status = status;
  }

  startTitleEdit(): void {
    this.titleEl.contentEditable = 'true';
    this.titleEl.focus();
    document.getSelection()?.selectAllChildren(this.titleEl);
  }

  startDescEdit(): void {
    this.descEl.dataset.editing = 'true';
    this.descEl.style.display = '';
    this.descEl.contentEditable = 'true';
    this.descEl.focus();
  }

  /** Dot centers sit at `GROUP_DOT_TOP + i * GROUP_DOT_STEP` — the same math as the editor's wire anchors. */
  renderDots(inCount: number, outCount: number): void {
    const fill = (host: HTMLElement, count: number): void => {
      host.textContent = '';
      for (let i = 0; i < count; i++) {
        const dot = document.createElement('span');
        dot.className = 'ff-group-dot';
        dot.style.top = `${GROUP_DOT_TOP + i * GROUP_DOT_STEP - 4.5}px`;
        host.appendChild(dot);
      }
    };
    fill(this.inDots, inCount);
    fill(this.outDots, outCount);
    const rows = Math.max(inCount, outCount);
    this.element.style.minHeight = rows > 0 ?
      `${GROUP_DOT_TOP + (rows - 1) * GROUP_DOT_STEP + 12}px` : '';
  }

  toDoc(): GroupDoc {
    return {
      id: this.id,
      title: this.title,
      description: this.description,
      memberIds: Array.from(this.memberIds),
      minimized: this.minimized,
      pos: {...this.pos},
    };
  }
}
