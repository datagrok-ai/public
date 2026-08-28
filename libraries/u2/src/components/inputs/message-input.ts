/* Compose box: a contenteditable editor whose value is Datagrok markup text — plain text carrying
   atomic, non-editable chips that serialize back to the exact platform token they were built from.
   Mentions come from `MentionProvider`s: the trigger wraps what was typed into a
   `span.u2-msg-query` and hangs the shared SuggestionList off it, so the popup follows the caret
   without caret-rect math; a pick swaps that span for a chip, and a restored value re-chips through
   the provider's token pattern. Platform-free — the users provider is src/dg/inputs/message-input.ts. */
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {AsyncSource, AsyncFetch, AsyncState} from '../../core/async-source.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import {ObjectRenderer} from '../../core/object-renderer.js';
import {button} from '../../core/elements.js';
import {iconButton} from '../actions/buttons.js';

export interface MentionProvider<T = unknown> {
  /** '@' this wave; the field makes '#' additive later. */
  trigger: string;
  fetch: AsyncFetch<T>;
  /** Rows through `listItem`/`markup`, chips through `markup`. */
  renderer?: ObjectRenderer<T>;
  /** Popup row; wins over the renderer's members. */
  renderItem?: (item: T) => HTMLElement;
  /** Chip text fallback. */
  caption(item: T): string;
  /** The token serialized into the value. */
  toMarkup(item: T): string;
  /** Regexp SOURCE matching this provider's tokens in a restored value; its capture groups reach
   * {@link renderToken}. Every provider's pattern is combined into one alternation, so it must
   * carry no backreferences and no group name another provider could also use. */
  tokenPattern?: string;
  renderToken?: (match: RegExpExecArray) => HTMLElement;
  minChars?: number;
  debounceMs?: number;
}

export interface MessageInputOptions {
  placeholder?: string;
  mentions?: MentionProvider<any>[];
  /** Default 'ctrlEnter' — the platform chat's own binding. */
  sendOn?: 'enter' | 'ctrlEnter';
  onSend?: (markup: string) => void | Promise<void>;
  /** Persists the unsent text under `u2-msg-draft-<key>`. */
  draftKey?: string;
  maxHeight?: number;
  /** Initial markup; wins over a stored draft. */
  value?: string;
}

interface Session {
  provider: MentionProvider<any>;
  span: HTMLElement;
  scope: Scope;
  source: AsyncSource<any>;
  list: SuggestionList<any>;
  query: Signal<string>;
  items: ReadonlySignal<any[]>;
}

const TRIGGER_ICONS: Record<string, string> = {'@': 'at', '#': 'hashtag'};
const BLOCK_TAGS = new Set(['DIV', 'P']);

export class MessageInput extends Control {
  readonly empty: ReadonlySignal<boolean>;
  readonly isMentionOpen: ReadonlySignal<boolean>;

  private readonly _options: MessageInputOptions;
  private readonly _providers: MentionProvider<any>[];
  private readonly _editor = document.createElement('div');
  private readonly _empty = signal(true);
  private readonly _open = signal(false);
  private readonly _draftKey: string | undefined;
  private readonly _tokens: RegExp | undefined;
  private readonly _owners: {provider: MentionProvider<any>, group: number}[] = [];
  private _session: Session | undefined;

  constructor(options: MessageInputOptions = {}) {
    super();
    this._options = options;
    this._providers = options.mentions ?? [];
    this._draftKey = options.draftKey === undefined ? undefined : `u2-msg-draft-${options.draftKey}`;
    this.empty = this._empty;
    this.isMentionOpen = this._open;
    this._tokens = this._buildTokens();

    this.root.classList.add('u2-msg');
    this.root.dataset.u2 = 'message-input';

    const editor = this._editor;
    editor.className = 'u2-msg-editor';
    editor.setAttribute('contenteditable', 'true');
    editor.setAttribute('role', 'textbox');
    editor.setAttribute('aria-multiline', 'true');
    editor.dataset.placeholder = options.placeholder ?? '';
    editor.style.maxHeight = `${options.maxHeight ?? 160}px`;

    this.run(() => {
      const tools = document.createElement('div');
      tools.className = 'u2-msg-tools';
      for (const provider of this._providers) {
        const trigger = iconButton(TRIGGER_ICONS[provider.trigger] ?? 'at',
          () => this.openMentions(provider), {tooltip: `Mention (${provider.trigger})`});
        // the default mousedown would blur the editor, and the session dies with the caret it needs
        this._listen(trigger, 'mousedown', (e) => e.preventDefault());
        tools.append(trigger);
      }
      const send = button('Send', () => this.send(), {primary: true});
      send.classList.add('u2-msg-send');
      const toolbar = document.createElement('div');
      toolbar.className = 'u2-msg-toolbar';
      toolbar.append(tools, send);
      this.root.append(editor, toolbar);
      this.effect(() => send.disabled = this._empty.value);
    });

    this._listen(editor, 'input', () => this._onInput());
    this._listen(editor, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(editor, 'blur', () => this._closeSession());
    this.own(() => this._closeSession());

    const initial = options.value ?? this._restoreDraft();
    if (initial !== undefined)
      this.value = initial;
    else
      this._refreshEmpty();
  }

  get value(): string {
    return MessageInput._serialize(this._editor).replace(/\n$/, '');
  }

  set value(markup: string) {
    this._endSession();
    this._editor.textContent = '';
    for (const node of this._parse(markup))
      this._editor.append(node);
    this._refreshEmpty();
  }

  insertMention(item: unknown, provider: MentionProvider<any> | undefined = this._providers[0]): void {
    if (provider === undefined)
      return;
    const session = this._session;
    if (session && session.provider === provider) {
      this._replaceQuery(session, item);
      return;
    }
    this._insertAtCaret(this._chipFor(provider, item), document.createTextNode(' '));
    this._afterEdit();
  }

  openMentions(provider: MentionProvider<any> | undefined = this._providers[0]): void {
    if (provider === undefined)
      return;
    this._closeSession();
    this._editor.focus();
    const span = document.createElement('span');
    span.className = 'u2-msg-query';
    span.textContent = provider.trigger;
    this._insertAtCaret(span);
    this._startSession(provider, span);
    this._caretTo(span.firstChild ?? span, provider.trigger.length);
    this._refreshEmpty();
  }

  /** The box is emptied only once the handler has accepted the message — a rejected send leaves the
   * text (and its draft) where the user can try again. */
  send(): void {
    if (this._empty.peek())
      return;
    const markup = this.value;
    const onSend = this._options.onSend;
    if (onSend === undefined) {
      this.clear();
      return;
    }
    void Promise.resolve(onSend(markup))
      .then(() => this.clear())
      .catch((e) => console.error('u2: message not sent', e));
  }

  clear(): void {
    this._closeSession();
    this._editor.textContent = '';
    this._afterEdit();
  }

  focus(): void {
    this._editor.focus();
  }

  // ---- serialization ------------------------------------------------------

  private static _serialize(node: Node): string {
    let out = '';
    for (const child of Array.from(node.childNodes)) {
      if (child.nodeType === 3) {
        out += child.textContent ?? '';
        continue;
      }
      if (child.nodeType !== 1)
        continue;
      const el = child as HTMLElement;
      if (el.classList.contains('u2-msg-chip'))
        out += el.dataset.token ?? '';
      else if (el.tagName === 'BR')
        out += '\n';
      else if (BLOCK_TAGS.has(el.tagName))
        out += `${MessageInput._serialize(el)}\n`;
      else
        out += MessageInput._serialize(el);
    }
    return out;
  }

  /** One alternation over every provider's token pattern; the wrapper group each alternative sits
   * in tells which provider a match belongs to. */
  private _buildTokens(): RegExp | undefined {
    const parts: string[] = [];
    let group = 1;
    for (const provider of this._providers) {
      const pattern = provider.tokenPattern;
      if (pattern === undefined)
        continue;
      parts.push(`(${pattern})`);
      this._owners.push({provider, group});
      group += MessageInput._groupCount(pattern) + 1;
    }
    return parts.length > 0 ? new RegExp(parts.join('|'), 'g') : undefined;
  }

  private static _groupCount(source: string): number {
    return new RegExp(`${source}|`).exec('')!.length - 1;
  }

  private _parse(markup: string): Node[] {
    const out: Node[] = [];
    const tokens = this._tokens;
    if (tokens === undefined) {
      MessageInput._pushText(out, markup);
      return out;
    }
    tokens.lastIndex = 0;
    let at = 0;
    for (let m = tokens.exec(markup); m !== null; m = tokens.exec(markup)) {
      MessageInput._pushText(out, markup.substring(at, m.index));
      out.push(this._restoreChip(m));
      at = m.index + m[0].length;
    }
    MessageInput._pushText(out, markup.substring(at));
    return out;
  }

  private static _pushText(out: Node[], text: string): void {
    if (text === '')
      return;
    const lines = text.split('\n');
    for (let i = 0; i < lines.length; i++) {
      if (i > 0)
        out.push(document.createElement('br'));
      if (lines[i] !== '')
        out.push(document.createTextNode(lines[i]));
    }
  }

  private _restoreChip(match: RegExpExecArray): HTMLElement {
    const owner = this._owners.find((o) => match[o.group] !== undefined);
    const provider = owner?.provider;
    const own = provider?.tokenPattern === undefined ? null :
      new RegExp(provider.tokenPattern).exec(match[0]);
    const content = provider?.renderToken && own ? provider.renderToken(own) :
      document.createTextNode(MessageInput._tokenText(own ?? match));
    return MessageInput._chip(match[0], content);
  }

  private static _tokenText(match: RegExpExecArray): string {
    const captured = match.slice(1).filter((g) => g !== undefined);
    return captured.length > 0 ? captured[captured.length - 1] : match[0];
  }

  private _chipFor(provider: MentionProvider<any>, item: unknown): HTMLElement {
    const markup = provider.renderer?.markup;
    const content = markup ? markup.call(provider.renderer, item) :
      document.createTextNode(provider.caption(item));
    return MessageInput._chip(provider.toMarkup(item), content);
  }

  private static _chip(token: string, content: Node): HTMLElement {
    const chip = document.createElement('span');
    chip.className = 'u2-msg-chip';
    chip.setAttribute('contenteditable', 'false');
    chip.dataset.u2 = 'msg-chip';
    chip.dataset.token = token;
    chip.append(content);
    return chip;
  }

  // ---- caret --------------------------------------------------------------

  /** The one caret read: a collapsed position inside the editor, or null where the host has no
   * selection to offer (headless) or it sits elsewhere. */
  private _selection(): {node: Node, offset: number} | null {
    if (typeof window.getSelection !== 'function')
      return null;
    const selection = window.getSelection();
    if (!selection || selection.rangeCount === 0)
      return null;
    const node = selection.anchorNode;
    if (!node || !this._editor.contains(node))
      return null;
    return {node, offset: selection.anchorOffset};
  }

  /** The one caret write; a host without Range (headless) simply keeps its own position. */
  private _caretTo(node: Node, offset: number): void {
    if (typeof window.getSelection !== 'function' || typeof document.createRange !== 'function')
      return;
    const selection = window.getSelection();
    if (!selection)
      return;
    const range = document.createRange();
    range.setStart(node, offset);
    range.collapse(true);
    selection.removeAllRanges();
    selection.addRange(range);
  }

  /** Inserts at the caret, splitting the text node it sits in; with no selection the nodes land at
   * the end of the editor. The caret ends up after the last one. */
  private _insertAtCaret(...nodes: Node[]): void {
    const at = this._selection();
    let parent: Node = this._editor;
    let before: Node | null = null;
    if (at && at.node.nodeType === 3) {
      const text = at.node;
      const data = text.textContent ?? '';
      parent = text.parentNode ?? this._editor;
      text.textContent = data.substring(0, at.offset);
      before = text.nextSibling;
      const tail = data.substring(at.offset);
      if (tail !== '')
        before = parent.insertBefore(document.createTextNode(tail), before);
    } else if (at) {
      parent = at.node;
      before = parent.childNodes[at.offset] ?? null;
    }
    for (const node of nodes)
      parent.insertBefore(node, before);
    const last = nodes[nodes.length - 1];
    if (last)
      this._caretTo(parent, Array.from(parent.childNodes).indexOf(last as ChildNode) + 1);
  }

  // ---- mention sessions ---------------------------------------------------

  private _startSession(provider: MentionProvider<any>, span: HTMLElement): void {
    const scope = new Scope();
    const query = signal(MessageInput._queryText(span, provider.trigger));
    const source = new AsyncSource<any>(provider.fetch, {debounceMs: provider.debounceMs});
    scope.own(() => source.dispose());
    const view = source.state;
    const items = computed<any[]>(() => {
      const state: AsyncState<any> = view.value;
      return state.kind === 'ready' ? state.items : [];
    });
    const renderer = provider.renderer;
    const listItem = renderer?.listItem?.bind(renderer) ?? renderer?.markup?.bind(renderer);
    const render = provider.renderItem ?? listItem ??
      ((item: any) => SuggestionList.row('u2-msg-text', provider.caption(item)));
    const list = new SuggestionList<any>({
      prefix: 'u2-msg',
      anchor: span,
      scope,
      items,
      view,
      text: query,
      minChars: provider.minChars ?? 0,
      render,
      autoHighlight: true,
      onPick: (index) => this._pick(index),
      onDismiss: () => this._closeSession(),
      onRetry: () => source.retry(),
    });
    const session = {provider, span, scope, source, list, query, items};
    this._session = session;
    this._open.value = true;
    list.open();
    this._refreshQuery(session);
  }

  private static _queryText(span: HTMLElement, trigger: string): string {
    return (span.textContent ?? '').substring(trigger.length);
  }

  private _refreshQuery(session: Session): void {
    const query = MessageInput._queryText(session.span, session.provider.trigger);
    session.query.value = query;
    if (query.length >= (session.provider.minChars ?? 0))
      session.source.query(query);
  }

  private _endSession(): Session | undefined {
    const session = this._session;
    if (session === undefined)
      return undefined;
    this._session = undefined;
    this._open.value = false;
    session.list.dismiss();
    session.scope.dispose();
    return session;
  }

  /** Dismissal: the query span goes back to the plain text it wrapped. */
  private _closeSession(): void {
    const session = this._endSession();
    const span = session?.span;
    if (!span)
      return;
    if (span.isConnected) {
      const parent = span.parentNode!;
      const text = document.createTextNode(span.textContent ?? '');
      parent.insertBefore(text, span);
      span.remove();
      this._caretTo(text, (text.textContent ?? '').length);
    }
    this._refreshEmpty();
  }

  private _pick(index: number): void {
    const session = this._session;
    if (session === undefined)
      return;
    const item = session.items.peek()[index];
    if (item === undefined)
      return;
    this._replaceQuery(session, item);
  }

  private _replaceQuery(session: Session, item: unknown): void {
    const {provider, span} = session;
    // ends the session first: the overlay watches its anchor and would fire a dismissal at us
    // in the middle of the swap
    this._endSession();
    const parent = span.parentNode;
    if (parent === null)
      return;
    const space = document.createTextNode(' ');
    parent.insertBefore(this._chipFor(provider, item), span);
    parent.insertBefore(space, span);
    span.remove();
    this._caretTo(space, 1);
    this._editor.focus();
    this._afterEdit();
  }

  /** Opens a session when the text before the caret ends in a provider's trigger word. */
  private _maybeOpen(): void {
    const at = this._selection();
    if (!at || at.node.nodeType !== 3)
      return;
    const before = (at.node.textContent ?? '').substring(0, at.offset);
    for (const provider of this._providers) {
      const match = new RegExp(`(?:^|\\s)(${MessageInput._escape(provider.trigger)}\\w*)$`).exec(before);
      if (match === null)
        continue;
      this._wrapQuery(provider, at.node, at.offset - match[1].length, at.offset);
      return;
    }
  }

  private static _escape(text: string): string {
    return text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  }

  private _wrapQuery(provider: MentionProvider<any>, text: Node, start: number, end: number): void {
    const parent = text.parentNode;
    if (parent === null)
      return;
    const data = text.textContent ?? '';
    const span = document.createElement('span');
    span.className = 'u2-msg-query';
    span.textContent = data.substring(start, end);
    text.textContent = data.substring(0, start);
    const next = text.nextSibling;
    parent.insertBefore(span, next);
    const tail = data.substring(end);
    if (tail !== '')
      parent.insertBefore(document.createTextNode(tail), next);
    this._startSession(provider, span);
    this._caretTo(span.firstChild ?? span, (span.textContent ?? '').length);
  }

  // ---- editing ------------------------------------------------------------

  private _onInput(): void {
    const session = this._session;
    if (session)
      this._refreshQuery(session);
    else
      this._maybeOpen();
    this._afterEdit();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const session = this._session;
    if (session) {
      switch (e.key) {
        case 'ArrowDown':
          e.preventDefault();
          session.list.move(1);
          return;
        case 'ArrowUp':
          e.preventDefault();
          session.list.move(-1);
          return;
        case 'Escape':
          e.stopPropagation();
          this._closeSession();
          return;
        case 'Enter':
          if (session.list.activeIndex.peek() >= 0) {
            e.preventDefault();
            e.stopPropagation();
            this._pick(session.list.activeIndex.peek());
            return;
          }
          break;
      }
    }
    switch (e.key) {
      case 'Enter':
        this._onEnter(e);
        break;
      case 'Backspace':
        this._onDelete(e, -1);
        break;
      case 'Delete':
        this._onDelete(e, 1);
        break;
    }
  }

  private _onEnter(e: KeyboardEvent): void {
    e.preventDefault();
    const modified = e.ctrlKey || e.metaKey;
    const sends = (this._options.sendOn ?? 'ctrlEnter') === 'enter' ? !modified && !e.shiftKey : modified;
    if (sends) {
      this.send();
      return;
    }
    this._insertAtCaret(document.createTextNode('\n'));
    this._afterEdit();
  }

  private _onDelete(e: KeyboardEvent, direction: -1 | 1): void {
    const chip = this._chipAtCaret(direction);
    if (chip === null)
      return;
    e.preventDefault();
    chip.remove();
    this._afterEdit();
  }

  private _chipAtCaret(direction: -1 | 1): HTMLElement | null {
    const at = this._selection();
    if (at === null)
      return null;
    let node: Node | null;
    if (at.node.nodeType === 3) {
      const length = (at.node.textContent ?? '').length;
      if (at.offset !== (direction < 0 ? 0 : length))
        return null;
      node = direction < 0 ? MessageInput._previous(at.node) : at.node.nextSibling;
    } else
      node = at.node.childNodes[direction < 0 ? at.offset - 1 : at.offset] ?? null;
    return node !== null && node.nodeType === 1 &&
      (node as HTMLElement).classList.contains('u2-msg-chip') ? node as HTMLElement : null;
  }

  private static _previous(node: Node): Node | null {
    const parent = node.parentNode;
    if (parent === null)
      return null;
    const index = Array.from(parent.childNodes).indexOf(node as ChildNode);
    return index > 0 ? parent.childNodes[index - 1] : null;
  }

  private _afterEdit(): void {
    const markup = this.value;
    this._refreshEmpty(markup);
    this._saveDraft(markup);
  }

  private _refreshEmpty(markup: string = this.value): void {
    this._empty.value = markup.trim() === '';
    this._editor.classList.toggle('u2-msg-editor-empty', markup === '');
  }

  // ---- draft --------------------------------------------------------------

  /** localStorage is unavailable in some hosts (headless tests, sandboxed frames) — then a draft
   * simply never survives. */
  private _saveDraft(markup: string): void {
    const key = this._draftKey;
    if (key === undefined)
      return;
    try {
      if (markup === '')
        window.localStorage.removeItem(key);
      else
        window.localStorage.setItem(key, markup);
    } catch {
      return;
    }
  }

  private _restoreDraft(): string | undefined {
    const key = this._draftKey;
    if (key === undefined)
      return undefined;
    try {
      return window.localStorage.getItem(key) ?? undefined;
    } catch {
      return undefined;
    }
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
