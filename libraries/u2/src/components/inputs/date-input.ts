/* Hand-written calendar machine (PLAN.md D6; Zag rejected), APG date-picker-dialog pattern.
   INVARIANT: `DateField` below is the single calendar implementation — `DateInput` and
   `DateTimeInput` are thin subclasses that differ only in `withTime`, so the field, the popup,
   the grid and the keyboard model are never forked. */
import {signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {button} from '../../core/elements.js';
import {iconButton} from '../actions/buttons.js';

export interface DateInputOptions extends InputOptions<Date | null> {
  /** Inclusive: days outside are disabled in the calendar, typed values clamp on commit. */
  min?: Date;
  max?: Date;
  /** 0 Sunday, 1 Monday (the default). */
  firstDayOfWeek?: 0 | 1;
}

export type DateState = 'idle' | 'focused' | 'open';

type DateEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'open'} | {type: 'toggleOpen'} |
  {type: 'escape'} | {type: 'dismiss'} | {type: 'today'} | {type: 'select', date: Date};

const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
  'September', 'October', 'November', 'December'];
const WEEKDAYS = ['Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa'];
const DATE = /^(\d{4})-(\d{1,2})-(\d{1,2})$/;
const DATE_TIME = /^(\d{4})-(\d{1,2})-(\d{1,2})[ T](\d{1,2}):(\d{1,2})$/;

function pad2(value: number): string {
  return value < 10 ? `0${value}` : String(value);
}

function startOfDay(date: Date): Date {
  return new Date(date.getFullYear(), date.getMonth(), date.getDate());
}

function startOfMonth(date: Date): Date {
  return new Date(date.getFullYear(), date.getMonth(), 1);
}

function addDays(date: Date, days: number): Date {
  return new Date(date.getFullYear(), date.getMonth(), date.getDate() + days);
}

function addMonths(date: Date, months: number): Date {
  const target = new Date(date.getFullYear(), date.getMonth() + months, 1);
  const last = new Date(target.getFullYear(), target.getMonth() + 1, 0).getDate();
  return new Date(target.getFullYear(), target.getMonth(), Math.min(date.getDate(), last));
}

function isoDay(date: Date): string {
  return `${date.getFullYear()}-${pad2(date.getMonth() + 1)}-${pad2(date.getDate())}`;
}

function fromIsoDay(iso: string): Date {
  const parts = iso.split('-').map(Number);
  return new Date(parts[0], parts[1] - 1, parts[2]);
}

abstract class DateField extends Input<Date | null, DateInputOptions> {
  private static _seq = 0;

  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<DateState>;

  // every field below is built by createEditor(), which the base constructor calls before
  // subclass initializers would run — none of them may carry an inline initializer
  private _id!: string;
  private _input!: HTMLInputElement;
  private _field!: HTMLElement;
  private _toggle!: HTMLElement;
  private _popup!: HTMLElement;
  private _title!: HTMLElement;
  private _weeks!: HTMLElement;
  private _hours: HTMLInputElement | undefined;
  private _minutes: HTMLInputElement | undefined;
  private _bad!: Signal<boolean>;
  private _open!: Signal<boolean>;
  private _machine!: Signal<DateState>;
  private _month!: Signal<Date>;
  private _active!: Signal<Date>;
  private _time!: Signal<number>;
  private _dayEls!: Signal<HTMLElement[]>;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: DateInputOptions) {
    super(options, null);
    this.isOpen = this._open;
    this.machineState = this._machine;
  }

  protected abstract get withTime(): boolean;

  protected createEditor(): HTMLElement {
    this._id = `u2-date-${++DateField._seq}`;
    this._bad = signal(false);
    this._open = signal(false);
    this._machine = signal<DateState>('idle');
    this._month = signal(startOfMonth(new Date()));
    this._active = signal(startOfDay(new Date()));
    this._time = signal(0);
    this._dayEls = signal<HTMLElement[]>([]);

    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = this._pattern;
    input.setAttribute('aria-label', this.options.label ?? this._pattern);

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (!echo)
        input.value = this._format(value);
    });
    this._listen(input, 'input', () => {
      const parsed = this._parse(input.value);
      this._bad.value = parsed === undefined;
      if (parsed === undefined)
        return;
      echo = true;
      try {
        this.value.value = parsed;
      } finally {
        echo = false;
      }
    });
    this._listen(input, 'blur', () => {
      this._commit();
      this._send({type: 'blur'});
    });
    this._listen(input, 'focus', () => this._send({type: 'focus'}));
    this._listen(input, 'keydown', (e) => this._onFieldKeyDown(e as KeyboardEvent));

    // `calendar-alt` on the options rail, always visible — the platform's own affordance
    // (`date_input.dart:19-28`, rail metrics `ui.css:239-242`)
    this._toggle = iconButton('calendar-alt', () => this._send({type: 'toggleOpen'}),
      {tooltip: 'Pick a date'});
    this._toggle.classList.add('u2-date-toggle');
    this._toggle.setAttribute('aria-haspopup', 'dialog');
    this.addOptions(this._toggle);

    const field = document.createElement('div');
    this._field = field;
    field.className = 'u2-date';
    field.append(input);
    this._popup = this._buildPopup();

    this.effect(() => {
      const open = this._open.value;
      this._toggle.setAttribute('aria-expanded', String(open));
      field.classList.toggle('u2-date-open', open);
    });
    this.effect(() => this._renderCalendar());
    this.effect(() => this._applyDays());
    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? `Not a ${this.withTime ? 'date-time' : 'date'}` : null);
    return field;
  }

  private get _pattern(): string {
    return this.withTime ? 'yyyy-MM-dd HH:mm' : 'yyyy-MM-dd';
  }

  private get _firstDay(): number {
    return this.options.firstDayOfWeek ?? 1;
  }

  private _buildPopup(): HTMLElement {
    const popup = DateField._el('div', 'u2-date-popup');
    popup.setAttribute('role', 'dialog');
    popup.setAttribute('aria-label', this.withTime ? 'Choose date and time' : 'Choose date');

    this._title = DateField._el('div', 'u2-date-title');
    this._title.id = `${this._id}-title`;
    const header = DateField._el('div', 'u2-date-header');
    header.append(
      iconButton('chevron-left', () => this._shiftMonth(-1), {tooltip: 'Previous month'}),
      this._title,
      iconButton('chevron-right', () => this._shiftMonth(1), {tooltip: 'Next month'}),
      button('Today', () => this._send({type: 'today'})));

    const grid = DateField._el('div', 'u2-date-grid');
    grid.setAttribute('role', 'grid');
    grid.setAttribute('aria-labelledby', this._title.id);
    const weekdays = DateField._el('div', 'u2-date-weekdays');
    weekdays.setAttribute('role', 'row');
    for (let i = 0; i < 7; i++) {
      const cell = DateField._el('div', 'u2-date-weekday', WEEKDAYS[(this._firstDay + i) % 7]);
      cell.setAttribute('role', 'columnheader');
      weekdays.append(cell);
    }
    this._weeks = DateField._el('div', 'u2-date-weeks');
    this._weeks.setAttribute('role', 'rowgroup');
    grid.append(weekdays, this._weeks);

    popup.append(header, grid);
    if (this.withTime)
      popup.append(this._buildTime());

    this._listen(grid, 'keydown', (e) => this._onGridKeyDown(e as KeyboardEvent));
    this._listen(grid, 'click', (e) => this._onGridClick(e));
    this._listen(popup, 'keydown', (e) => this._onPopupKeyDown(e as KeyboardEvent));
    this._listen(popup, OVERLAY_CLOSE_EVENT, () => this._send({type: 'dismiss'}));
    return popup;
  }

  private _buildTime(): HTMLElement {
    const row = DateField._el('div', 'u2-date-time');
    this._hours = this._segment('Hours');
    this._minutes = this._segment('Minutes');
    row.append(DateField._el('span', 'u2-date-time-label', 'Time'), this._hours,
      DateField._el('span', 'u2-date-time-sep', ':'), this._minutes);
    return row;
  }

  private _segment(label: string): HTMLInputElement {
    const input = document.createElement('input');
    input.className = 'u2-date-time-seg';
    input.type = 'text';
    input.inputMode = 'numeric';
    input.autocomplete = 'off';
    input.maxLength = 2;
    input.setAttribute('aria-label', label);
    this._listen(input, 'input', () => this._onTimeInput(input));
    this._listen(input, 'keydown', (e) => this._onTimeKeyDown(e as KeyboardEvent, input));
    this._listen(input, 'blur', () => this._syncTime());
    return input;
  }

  private _send(ev: DateEvent): void {
    const state = this._machine.peek();
    switch (ev.type) {
      case 'focus':
        if (state === 'idle')
          this._machine.value = 'focused';
        break;
      case 'blur':
        if (state !== 'open')
          this._machine.value = 'idle';
        break;
      case 'open':
        this._openPopup();
        break;
      case 'toggleOpen':
        if (state === 'open')
          this._dismiss(true);
        else
          this._openPopup();
        break;
      case 'escape':
        this._dismiss(true);
        break;
      case 'dismiss':
        this._dismiss(false);
        break;
      case 'today':
        this._today();
        break;
      case 'select':
        this._select(ev.date);
        break;
    }
  }

  private _openPopup(): void {
    if (this._open.peek())
      return;
    const value = this.value.peek();
    const base = value ?? this._defaultDay();
    this._month.value = startOfMonth(base);
    this._active.value = startOfDay(base);
    this._time.value = value === null ? 0 : value.getHours() * 60 + value.getMinutes();
    this._machine.value = 'open';
    this._closeOverlay = Overlay.show(this._field, this._popup, this.scope);
    this._open.value = true;
    this._syncTime();
    this._focusActive();
  }

  private _dismiss(refocus: boolean): void {
    if (!this._open.peek())
      return;
    this._open.value = false;
    this._machine.value = refocus ? 'focused' : 'idle';
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
    if (refocus)
      this._input.focus();
  }

  private _defaultDay(): Date {
    const today = startOfDay(new Date());
    const {min, max} = this.options;
    if (min !== undefined && today < startOfDay(min))
      return startOfDay(min);
    if (max !== undefined && today > startOfDay(max))
      return startOfDay(max);
    return today;
  }

  private _shiftMonth(delta: number): void {
    this._goto(addMonths(this._active.peek(), delta), false);
  }

  private _today(): void {
    const today = startOfDay(new Date());
    if (this._disabled(today))
      this._goto(today, true);
    else
      this._select(today);
  }

  private _goto(date: Date, focus: boolean): void {
    if (date.getFullYear() !== this._month.peek().getFullYear() ||
        date.getMonth() !== this._month.peek().getMonth())
      this._month.value = startOfMonth(date);
    this._active.value = date;
    if (focus)
      this._focusActive();
  }

  private _select(date: Date): void {
    if (this._disabled(date))
      return;
    this._set(this._combine(date));
    this._dismiss(true);
  }

  private _combine(day: Date): Date {
    if (!this.withTime)
      return startOfDay(day);
    const time = this._time.peek();
    return new Date(day.getFullYear(), day.getMonth(), day.getDate(),
      Math.floor(time / 60), time % 60);
  }

  private _focusActive(): void {
    const iso = isoDay(this._active.peek());
    for (const el of this._dayEls.peek()) {
      if (el.dataset.date === iso)
        el.focus();
    }
  }

  private _onFieldKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Enter') {
      this._commit();
      return;
    }
    if (e.key !== 'ArrowDown')
      return;
    e.preventDefault();
    this._send({type: 'open'});
  }

  private _onPopupKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Escape')
      this._send({type: 'escape'});
  }

  private _onGridKeyDown(e: KeyboardEvent): void {
    const active = this._active.peek();
    let next: Date;
    switch (e.key) {
      case 'ArrowLeft':
        next = addDays(active, -1);
        break;
      case 'ArrowRight':
        next = addDays(active, 1);
        break;
      case 'ArrowUp':
        next = addDays(active, -7);
        break;
      case 'ArrowDown':
        next = addDays(active, 7);
        break;
      case 'PageUp':
        next = addMonths(active, e.shiftKey ? -12 : -1);
        break;
      case 'PageDown':
        next = addMonths(active, e.shiftKey ? 12 : 1);
        break;
      case 'Home':
        next = addDays(active, -this._weekOffset(active));
        break;
      case 'End':
        next = addDays(active, 6 - this._weekOffset(active));
        break;
      case 'Enter':
      case ' ':
        e.preventDefault();
        this._send({type: 'select', date: active});
        return;
      default:
        return;
    }
    e.preventDefault();
    this._goto(next, true);
  }

  private _onGridClick(e: Event): void {
    const cell = (e.target as HTMLElement).closest('.u2-date-day') as HTMLElement | null;
    if (cell)
      this._send({type: 'select', date: fromIsoDay(cell.dataset.date as string)});
  }

  private _onTimeInput(input: HTMLInputElement): void {
    const digits = input.value.replace(/\D/g, '').slice(0, 2);
    if (digits !== input.value)
      input.value = digits;
    this._readTime();
    if (digits.length === 2 && input === this._hours)
      this._minutes!.focus();
  }

  private _onTimeKeyDown(e: KeyboardEvent, input: HTMLInputElement): void {
    if (e.key === 'Enter') {
      e.preventDefault();
      this._send({type: 'select', date: this._active.peek()});
      return;
    }
    const step = e.key === 'ArrowUp' ? 1 : e.key === 'ArrowDown' ? -1 : 0;
    if (step === 0)
      return;
    e.preventDefault();
    const wrap = input === this._hours ? 24 : 60;
    input.value = pad2(((Number(input.value) || 0) + step + wrap) % wrap);
    this._readTime();
  }

  private _readTime(): void {
    const hours = Math.min(Number(this._hours!.value) || 0, 23);
    const minutes = Math.min(Number(this._minutes!.value) || 0, 59);
    this._time.value = hours * 60 + minutes;
    const value = this.value.peek();
    if (value !== null)
      this._set(this._combine(value));
  }

  private _syncTime(): void {
    if (!this._hours || !this._minutes)
      return;
    const time = this._time.peek();
    this._hours.value = pad2(Math.floor(time / 60));
    this._minutes.value = pad2(time % 60);
  }

  private _weekOffset(date: Date): number {
    return (date.getDay() - this._firstDay + 7) % 7;
  }

  private _renderCalendar(): void {
    if (!this._open.value) {
      this._dayEls.value = [];
      return;
    }
    const month = this._month.value;
    this._title.textContent = `${MONTHS[month.getMonth()]} ${month.getFullYear()}`;
    const first = addDays(month, -this._weekOffset(month));
    const today = isoDay(new Date());
    const els: HTMLElement[] = [];
    this._weeks.textContent = '';
    for (let week = 0; week < 6; week++) {
      const row = DateField._el('div', 'u2-date-week');
      row.setAttribute('role', 'row');
      for (let day = 0; day < 7; day++) {
        const date = addDays(first, week * 7 + day);
        const cell = this._dayCell(date, month, today);
        row.append(cell);
        els.push(cell);
      }
      this._weeks.append(row);
    }
    this._dayEls.value = els;
  }

  private _dayCell(date: Date, month: Date, today: string): HTMLElement {
    const cell = DateField._el('div', 'u2-date-day', String(date.getDate()));
    const iso = isoDay(date);
    cell.dataset.date = iso;
    cell.tabIndex = -1;
    cell.setAttribute('role', 'gridcell');
    cell.setAttribute('aria-label', `${MONTHS[date.getMonth()]} ${date.getDate()}, ${date.getFullYear()}`);
    if (date.getMonth() !== month.getMonth())
      cell.classList.add('u2-date-adjacent');
    if (iso === today) {
      cell.classList.add('u2-date-today');
      cell.setAttribute('aria-current', 'date');
    }
    if (this._disabled(date)) {
      cell.classList.add('u2-date-day-disabled');
      cell.setAttribute('aria-disabled', 'true');
    }
    return cell;
  }

  private _applyDays(): void {
    const value = this.value.value;
    const selected = value === null ? null : isoDay(value);
    const active = isoDay(this._active.value);
    for (const el of this._dayEls.value) {
      const on = el.dataset.date === selected;
      el.setAttribute('aria-selected', String(on));
      el.classList.toggle('u2-date-selected', on);
      el.tabIndex = el.dataset.date === active ? 0 : -1;
    }
  }

  private _commit(): void {
    const parsed = this._parse(this._input.value);
    if (parsed !== undefined)
      this._set(parsed);
  }

  private _set(value: Date | null): void {
    this._bad.value = false;
    this.value.value = value === null ? null : this._clamp(value);
    this._input.value = this._format(this.value.peek());
  }

  private _disabled(date: Date): boolean {
    const {min, max} = this.options;
    const day = startOfDay(date).getTime();
    return (min !== undefined && day < startOfDay(min).getTime()) ||
      (max !== undefined && day > startOfDay(max).getTime());
  }

  private _clamp(value: Date): Date {
    const {min, max} = this.options;
    if (min !== undefined && value < min)
      return this._normalize(min);
    if (max !== undefined && value > max)
      return this._normalize(max);
    return value;
  }

  private _normalize(date: Date): Date {
    if (!this.withTime)
      return startOfDay(date);
    return new Date(date.getFullYear(), date.getMonth(), date.getDate(),
      date.getHours(), date.getMinutes());
  }

  private _format(value: Date | null): string {
    if (value === null)
      return '';
    const day = isoDay(value);
    return this.withTime ? `${day} ${pad2(value.getHours())}:${pad2(value.getMinutes())}` : day;
  }

  private _parse(text: string): Date | null | undefined {
    const trimmed = text.trim();
    if (!trimmed)
      return null;
    const match = (this.withTime ? DATE_TIME : DATE).exec(trimmed);
    if (!match)
      return undefined;
    const month = Number(match[2]);
    const day = Number(match[3]);
    const hours = match[4] === undefined ? 0 : Number(match[4]);
    const minutes = match[5] === undefined ? 0 : Number(match[5]);
    if (month < 1 || month > 12 || day < 1 || hours > 23 || minutes > 59)
      return undefined;
    const date = new Date(Number(match[1]), month - 1, day, hours, minutes);
    return date.getMonth() === month - 1 && date.getDate() === day ? date : undefined;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }

  private static _el(tag: string, cls: string, text?: string): HTMLElement {
    const el = document.createElement(tag);
    el.className = cls;
    if (text !== undefined)
      el.textContent = text;
    return el;
  }
}

/** Date field in `yyyy-MM-dd`, committed at midnight local time. v1 is a single date with no
 * ranges and no locale switching — both are future work. */
export class DateInput extends DateField {
  constructor(options: DateInputOptions = {}) {
    super(options);
    this.root.dataset.u2 = 'date-input';
  }

  protected get withTime(): boolean {
    return false;
  }
}

/** Date-time field in `yyyy-MM-dd HH:mm` (24h, local), committed with seconds and milliseconds
 * zeroed. v1 is a single date-time with no ranges and no locale switching — both are future work. */
export class DateTimeInput extends DateField {
  constructor(options: DateInputOptions = {}) {
    super(options);
    this.root.dataset.u2 = 'datetime-input';
  }

  protected get withTime(): boolean {
    return true;
  }
}
