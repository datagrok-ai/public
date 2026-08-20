/* DateInput / DateTimeInput: the text contract (parse, format, validity, clamping), the popup
   machine, and the calendar grid — one calendar shared by both, so the grid tests run on either. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {DateInput, DateTimeInput} from '../src/components/inputs/date-input.js';

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(component) {
  document.body.append(component.root);
  return component;
}

function editor(input) {
  return input.root.querySelector('input');
}

function toggle(input) {
  return input.root.querySelector('.u2-date-toggle');
}

function popup() {
  return document.body.querySelector('.u2-date-popup');
}

function grid() {
  return popup().querySelector('.u2-date-grid');
}

function days() {
  return popup().querySelectorAll('.u2-date-day');
}

function day(iso) {
  return popup().querySelector(`.u2-date-day[data-date='${iso}']`);
}

function active() {
  return days().filter((d) => d.tabIndex === 0);
}

function type(input, text) {
  const el = editor(input);
  el.value = text;
  fire(el, 'input');
  return el;
}

async function open(input) {
  fire(toggle(input), 'click');
  await flush();
  return popup();
}

smoke('parse and format round-trip, single-digit month and day, empty is null', async () => {
  const input = mount(new DateInput({label: 'Date'}));
  assert.equal(input.value.value, null);
  assert.equal(editor(input).value, '');

  type(input, '2026-08-15');
  assert.equal(input.value.value.getFullYear(), 2026);
  assert.equal(input.value.value.getMonth(), 7);
  assert.equal(input.value.value.getDate(), 15);
  assert.equal(input.value.value.getHours(), 0, 'DateInput commits midnight local');
  assert.equal(input.validity.value, null);

  type(input, '2026-8-5');
  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '2026-08-05', 'commit reformats to the canonical pattern');
  assert.equal(input.value.value.getDate(), 5);

  type(input, '');
  assert.equal(input.value.value, null);
  assert.equal(input.validity.value, null);

  const dt = mount(new DateTimeInput({label: 'When'}));
  type(dt, '2026-08-15 7:5');
  fire(editor(dt), 'blur');
  assert.equal(editor(dt).value, '2026-08-15 07:05');
  assert.equal(dt.value.value.getHours(), 7);
  assert.equal(dt.value.value.getMinutes(), 5);
  assert.equal(dt.value.value.getSeconds(), 0);
  assert.equal(dt.value.value.getMilliseconds(), 0);
  input.dispose();
  dt.dispose();
});

smoke('unparseable text stays on screen, marks invalid and never writes the value', async () => {
  const input = mount(new DateInput({label: 'Date', value: new Date(2026, 7, 15)}));
  const good = input.value.value;

  type(input, '2026-13-40');
  assert.equal(input.validity.value, 'Not a date');
  assert.equal(input.value.value, good, 'the value signal keeps the last good date');
  assert.equal(editor(input).classList.contains('u2-invalid'), false);
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Not a date');

  type(input, '2026-02-30');
  assert.equal(input.validity.value, 'Not a date', 'February 30 does not roll over into March');
  assert.equal(input.value.value, good);

  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '2026-02-30', 'a bad commit keeps the raw text');
  assert.equal(input.validity.value, 'Not a date');

  type(input, '2026-08-16');
  assert.equal(input.validity.value, null);
  assert.equal(input.value.value.getDate(), 16);

  const dt = mount(new DateTimeInput({label: 'When'}));
  type(dt, '2026-08-15 25:00');
  assert.equal(dt.validity.value, 'Not a date-time');
  assert.equal(dt.value.value, null);
  input.dispose();
  dt.dispose();
});

smoke('commit clamps to min/max; programmatic writes are echo-suppressed', async () => {
  const min = new Date(2026, 7, 10);
  const max = new Date(2026, 7, 20);
  const input = mount(new DateInput({label: 'Date', min, max}));

  type(input, '2026-08-01');
  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '2026-08-10', 'below min clamps up on commit');
  assert.equal(input.value.value.getTime(), min.getTime());

  type(input, '2026-08-25');
  assert.equal(input.value.value.getDate(), 25, 'mid-typing writes are not clamped');
  fire(editor(input), 'keydown', {key: 'Enter'});
  assert.equal(editor(input).value, '2026-08-20');
  assert.equal(input.value.value.getTime(), max.getTime());

  input.value.value = new Date(2026, 7, 14);
  assert.equal(editor(input).value, '2026-08-14', 'a programmatic write echoes into the text');
  input.value.value = null;
  assert.equal(editor(input).value, '');
  input.dispose();
});

smoke('machine: toggle and ArrowDown open, Esc closes and refocuses the editor', async () => {
  const input = mount(new DateInput({label: 'Date'}));
  assert.equal(input.machineState.value, 'idle');
  assert.equal(input.isOpen.value, false);
  assert.equal(toggle(input).getAttribute('aria-expanded'), 'false');
  assert.equal(toggle(input).getAttribute('aria-haspopup'), 'dialog');

  editor(input).focus();
  assert.equal(input.machineState.value, 'focused');

  await open(input);
  assert.equal(input.isOpen.value, true);
  assert.equal(input.machineState.value, 'open');
  assert.equal(toggle(input).getAttribute('aria-expanded'), 'true');
  assert.equal(popup().getAttribute('role'), 'dialog');
  assert.equal(grid().getAttribute('role'), 'grid');
  assert.equal(popup().querySelectorAll('.u2-date-week').length, 6);
  assert.equal(days().length, 42);

  fire(grid(), 'keydown', {key: 'Escape'});
  assert.equal(input.isOpen.value, false);
  assert.equal(input.machineState.value, 'focused');
  assert.equal(popup(), null, 'the overlay is gone');
  assert.equal(document.activeElement, editor(input), 'focus returns to the editor');

  fire(editor(input), 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(input.isOpen.value, true, 'ArrowDown reopens');

  fire(toggle(input), 'click');
  assert.equal(input.isOpen.value, false, 'the toggle closes an open popup');
  input.dispose();
});

smoke('the grid marks today, the selection, adjacent months and out-of-range days', async () => {
  const today = new Date();
  const input = mount(new DateInput({
    label: 'Date',
    value: new Date(2026, 7, 15),
    min: new Date(2026, 7, 10),
    max: new Date(2026, 7, 20),
    firstDayOfWeek: 1,
  }));
  await open(input);

  assert.equal(popup().querySelector('.u2-date-title').textContent, 'August 2026');
  assert.deepEqual(popup().querySelectorAll('.u2-date-weekday').map((d) => d.textContent),
    ['Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa', 'Su']);
  assert.equal(days()[0].dataset.date, '2026-07-27', 'the week starts on Monday');

  const picked = day('2026-08-15');
  assert.equal(picked.getAttribute('aria-selected'), 'true');
  assert.equal(picked.classList.contains('u2-date-selected'), true);
  assert.equal(day('2026-08-14').getAttribute('aria-selected'), 'false');
  assert.equal(day('2026-07-31').classList.contains('u2-date-adjacent'), true);
  assert.equal(day('2026-08-09').getAttribute('aria-disabled'), 'true');
  assert.equal(day('2026-08-21').getAttribute('aria-disabled'), 'true');
  assert.equal(day('2026-08-10').hasAttribute('aria-disabled'), false);

  const current = days().filter((d) => d.getAttribute('aria-current') === 'date');
  const inMonth = today.getFullYear() === 2026 && today.getMonth() === 7;
  assert.equal(current.length, inMonth ? 1 : 0);

  fire(day('2026-08-09'), 'click');
  assert.equal(input.value.value.getDate(), 15, 'a disabled day is not selectable');
  assert.equal(input.isOpen.value, true);

  fire(day('2026-08-18'), 'click');
  assert.equal(input.value.value.getDate(), 18, 'a click selects');
  assert.equal(input.isOpen.value, false, '…and closes');
  assert.equal(editor(input).value, '2026-08-18');

  await open(input);
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'August 2026',
    'the visible month follows the value on reopen');
  input.dispose();
});

smoke('grid keyboard: arrows, Home/End, PageUp/PageDown, Shift+Page, Enter selects', async () => {
  const input = mount(new DateInput({label: 'Date', value: new Date(2026, 7, 15)}));
  await open(input);
  assert.equal(active().length, 1);
  assert.equal(active()[0].dataset.date, '2026-08-15', 'opening makes the picked day active');
  assert.equal(document.activeElement, active()[0], 'the roving tabindex cell has focus');

  fire(grid(), 'keydown', {key: 'ArrowRight'});
  assert.equal(active()[0].dataset.date, '2026-08-16');
  fire(grid(), 'keydown', {key: 'ArrowDown'});
  assert.equal(active()[0].dataset.date, '2026-08-23');
  fire(grid(), 'keydown', {key: 'ArrowUp'});
  fire(grid(), 'keydown', {key: 'ArrowLeft'});
  assert.equal(active()[0].dataset.date, '2026-08-15');

  fire(grid(), 'keydown', {key: 'Home'});
  assert.equal(active()[0].dataset.date, '2026-08-10', 'Home goes to the week start');
  fire(grid(), 'keydown', {key: 'End'});
  assert.equal(active()[0].dataset.date, '2026-08-16', 'End goes to the week end');

  fire(grid(), 'keydown', {key: 'PageUp'});
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'July 2026');
  assert.equal(active()[0].dataset.date, '2026-07-16');
  fire(grid(), 'keydown', {key: 'PageDown'});
  fire(grid(), 'keydown', {key: 'PageDown'});
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'September 2026');
  fire(grid(), 'keydown', {key: 'PageUp', shiftKey: true});
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'September 2025');
  assert.equal(input.value.value.getFullYear(), 2026, 'navigation never writes the value');

  fire(grid(), 'keydown', {key: 'PageDown', shiftKey: true});
  fire(grid(), 'keydown', {key: 'Enter'});
  assert.equal(input.isOpen.value, false);
  assert.equal(editor(input).value, '2026-09-16');
  assert.equal(input.value.value.getMonth(), 8);
  input.dispose();
});

smoke('header: prev/next month move the calendar, Today picks today', async () => {
  const input = mount(new DateInput({label: 'Date', value: new Date(2026, 7, 15)}));
  await open(input);
  const [prev, next] = popup().querySelectorAll('.u2-date-header .u2-icon-btn');

  fire(prev, 'click');
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'July 2026');
  fire(next, 'click');
  fire(next, 'click');
  assert.equal(popup().querySelector('.u2-date-title').textContent, 'September 2026');
  assert.equal(input.value.value.getMonth(), 7, 'paging never writes the value');

  fire(popup().querySelector('.u2-date-header .u2-btn:not(.u2-icon-btn)'), 'click');
  const today = new Date();
  assert.equal(input.isOpen.value, false);
  assert.equal(input.value.value.getDate(), today.getDate());
  assert.equal(input.value.value.getMonth(), today.getMonth());
  assert.equal(input.value.value.getHours(), 0);
  input.dispose();
});

smoke('DateTimeInput: segments step and wrap, and a day change keeps the time', async () => {
  const input = mount(new DateTimeInput({label: 'When', value: new Date(2026, 7, 15, 9, 30)}));
  assert.equal(editor(input).value, '2026-08-15 09:30');
  await open(input);

  const [hours, minutes] = popup().querySelectorAll('.u2-date-time-seg');
  assert.equal(hours.value, '09', 'the segments seed from the value');
  assert.equal(minutes.value, '30');

  fire(hours, 'keydown', {key: 'ArrowUp'});
  assert.equal(hours.value, '10');
  assert.equal(input.value.value.getHours(), 10, 'stepping commits into the value');
  hours.value = '23';
  fire(hours, 'input');
  fire(hours, 'keydown', {key: 'ArrowUp'});
  assert.equal(hours.value, '00', 'hours wrap');
  minutes.value = '00';
  fire(minutes, 'input');
  fire(minutes, 'keydown', {key: 'ArrowDown'});
  assert.equal(minutes.value, '59', 'minutes wrap');
  assert.equal(editor(input).value, '2026-08-15 00:59');

  fire(day('2026-08-18'), 'click');
  assert.equal(editor(input).value, '2026-08-18 00:59', 'picking a day keeps the time of day');
  assert.equal(input.value.value.getMinutes(), 59);
  assert.equal(input.isOpen.value, false);

  await open(input);
  const segments = popup().querySelectorAll('.u2-date-time-seg');
  segments[0].value = '14';
  fire(segments[0], 'input');
  assert.equal(document.activeElement, segments[1], 'two digits auto-advance to the minutes');
  segments[1].value = '05';
  fire(segments[1], 'input');
  fire(segments[1], 'keydown', {key: 'Enter'});
  assert.equal(input.isOpen.value, false, 'Enter commits date and time together');
  assert.equal(editor(input).value, '2026-08-18 14:05');
  input.dispose();
});

smoke('firstDayOfWeek: 0 starts the week on Sunday', async () => {
  const input = mount(new DateInput({label: 'Date', value: new Date(2026, 7, 15), firstDayOfWeek: 0}));
  await open(input);
  assert.deepEqual(popup().querySelectorAll('.u2-date-weekday').map((d) => d.textContent),
    ['Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa']);
  assert.equal(days()[0].dataset.date, '2026-07-26');

  fire(grid(), 'keydown', {key: 'Home'});
  assert.equal(active()[0].dataset.date, '2026-08-09', 'the week starts on Sunday');
  input.dispose();
});

smoke('two inputs on one signal stay in step; dispose closes the popup and kills the listeners',
  async () => {
    const shared = mount(new DateInput({label: 'From', value: new Date(2026, 7, 15)}));
    const other = mount(new DateInput({label: 'Also from', bind: shared.value}));
    assert.equal(editor(other).value, '2026-08-15');
    type(shared, '2026-08-20');
    assert.equal(editor(other).value, '2026-08-20');
    other.dispose();

    const before = Scope.liveCount;
    const input = mount(new DateInput({label: 'Date'}));
    await open(input);
    const opened = popup();
    assert.ok(Scope.liveCount > before);

    input.dispose();
    assert.equal(opened.isConnected, false, 'dispose closes the popup');
    assert.equal(Scope.liveCount, before);

    fire(toggle(input), 'click');
    await flush();
    assert.equal(popup(), null, 'listeners died with the scope');
    shared.dispose();
  });
