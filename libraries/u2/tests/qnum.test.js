/* QNum: the packed-double contract ported from ddt/lib/src/utils/qnum.dart:17-93, pinned at byte
   level, plus the editor over it.

   Every fixture below is derived from the Dart source by hand, not captured from this port: the
   double's own bit pattern is IEEE-754 (identical in Dart and JS), and the packing is Dart's
   `byteData.setInt8(7, (last & 0xFC) | q)` on the big-endian byte 7 — the low mantissa byte, since
   Dart's ByteData and the DOM's DataView both default to big-endian. The derivation of each
   fixture's last byte is written out where it is used. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {QNum} from '../src/core/qnum.js';
import {QNumInput} from '../src/components/inputs/qnum-input.js';

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

/** The eight big-endian bytes of a double, read through a DataView of this test's own. */
function hex(x) {
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, x);
  return [...new Uint8Array(view.buffer)].map((b) => b.toString(16).padStart(2, '0')).join('');
}

function editor(input) {
  return input.root.querySelector('input');
}

function type(input, text) {
  const el = editor(input);
  el.value = text;
  fire(el, 'input');
}

test('the qualifier constants are Dart\'s (qnum.dart:18-20)', () => {
  assert.equal(QNum.LESS, 1);
  assert.equal(QNum.EXACT, 2);
  assert.equal(QNum.GREATER, 3);
});

test('fixture 5.2 — a positive double whose low byte already carries bits', () => {
  // 5.2 is 4014cccccccccccd; last byte cd = 1100 1101, so (cd & fc) = cc and the three
  // qualifiers land on cd (|1), ce (|2), cf (|3).
  assert.equal(hex(5.2), '4014cccccccccccd');
  assert.equal(hex(QNum.create(5.2, QNum.LESS)), '4014cccccccccccd');
  assert.equal(hex(QNum.create(5.2, QNum.EXACT)), '4014ccccccccccce');
  assert.equal(hex(QNum.create(5.2)), '4014ccccccccccce', 'EXACT is the default (qnum.dart:49)');
  assert.equal(hex(QNum.create(5.2, QNum.GREATER)), '4014cccccccccccf');

  for (const q of [QNum.LESS, QNum.EXACT, QNum.GREATER]) {
    const packed = QNum.create(5.2, q);
    assert.equal(QNum.getQ(packed), q);
    // getValue clears the two bits: cd/ce/cf all come back as cc
    assert.equal(hex(QNum.getValue(packed)), '4014cccccccccccc');
    assert.equal(QNum.getValue(packed), 5.199999999999999, 'the two borrowed bits are gone for good');
  }
});

test('fixture -12.75 — a negative double with a clear low byte, so packing is lossless', () => {
  // -12.75 is c029800000000000; last byte 00, so (00 & fc) | q is just q, and clearing gives
  // back the exact same double.
  assert.equal(hex(-12.75), 'c029800000000000');
  assert.equal(hex(QNum.create(-12.75, QNum.LESS)), 'c029800000000001');
  assert.equal(hex(QNum.create(-12.75, QNum.EXACT)), 'c029800000000002');
  assert.equal(hex(QNum.create(-12.75, QNum.GREATER)), 'c029800000000003');

  const packed = QNum.create(-12.75, QNum.GREATER);
  assert.equal(QNum.getQ(packed), QNum.GREATER);
  assert.equal(QNum.getValue(packed), -12.75, 'nothing to lose in a value with zero low bits');
});

test('fixture 5e-324 — the smallest subnormal, where the qualifier is the whole number', () => {
  // 5e-324 is 0000000000000001: the entire magnitude lives in the two bits the qualifier takes.
  // (01 & fc) | 2 = 02, which is 1e-323, and clearing leaves zero.
  assert.equal(hex(5e-324), '0000000000000001');
  assert.equal(hex(QNum.create(5e-324, QNum.EXACT)), '0000000000000002');
  assert.equal(QNum.getValue(QNum.create(5e-324, QNum.EXACT)), 0);

  // a subnormal one byte up keeps everything: 1e-320 is 00000000000007e8, and e8 & fc = e8
  assert.equal(hex(1e-320), '00000000000007e8');
  assert.equal(hex(QNum.create(1e-320, QNum.LESS)), '00000000000007e9');
  assert.equal(QNum.getValue(QNum.create(1e-320, QNum.LESS)), 1e-320);
  assert.equal(QNum.getQ(QNum.create(1e-320, QNum.LESS)), QNum.LESS);
});

test('zero: the qualifier is all that is left, and it reads back', () => {
  assert.equal(hex(QNum.create(0, QNum.LESS)), '0000000000000001');
  assert.equal(QNum.getQ(QNum.create(0, QNum.LESS)), QNum.LESS);
  assert.equal(QNum.getValue(QNum.create(0, QNum.GREATER)), 0);
  assert.equal(QNum.getQ(0), 0, 'a plain zero has no qualifier bits at all');
  assert.equal(QNum.qualifier(0), '', 'and reads as exact (qnum.dart:86)');
});

test('parse: the optional < / > prefix, and null for anything unparseable (qnum.dart:62-83)', () => {
  const less = QNum.parse('<5.2');
  assert.equal(QNum.getQ(less), QNum.LESS);
  assert.equal(QNum.getValue(less), 5.199999999999999);
  assert.equal(hex(less), '4014cccccccccccd');

  assert.equal(QNum.getQ(QNum.parse('>1e3')), QNum.GREATER);
  assert.equal(QNum.getValue(QNum.parse('>1e3')), 1000);
  assert.equal(QNum.getQ(QNum.parse('10')), QNum.EXACT, 'no prefix is exact');
  assert.equal(QNum.getQ(QNum.parse('  <  0.5')), QNum.LESS, 'leading spaces are skipped either side');
  assert.equal(QNum.getValue(QNum.parse('-12.75')), -12.75);

  assert.equal(QNum.parse(''), null);
  assert.equal(QNum.parse('   '), null);
  assert.equal(QNum.parse('<'), null);
  assert.equal(QNum.parse('abc'), null);
  assert.equal(QNum.parse('5,2'), null);
});

test('format: the glyph in front, and the packing rounded back off (qnum.dart:92)', () => {
  assert.equal(QNum.format(QNum.create(5.2, QNum.LESS)), '<5.2');
  assert.equal(QNum.format(QNum.create(5.2)), '5.2', 'exact shows no glyph');
  assert.equal(QNum.format(QNum.create(1000, QNum.GREATER)), '>1000');
  assert.equal(QNum.format(null), '');
  assert.equal(QNum.format(QNum.create(5.2, QNum.LESS), (v) => v.toFixed(1)), '<5.2',
    'a format takes over the number, never the glyph');
  assert.throws(() => QNum.qualifier(7), /Unknown qualifier type: 7/);
});

test('parse and format round-trip through each other', () => {
  for (const text of ['<5.2', '>1e3', '10', '-12.75', '<0.001'])
    assert.equal(QNum.format(QNum.parse(text)), text === '>1e3' ? '>1000' : text);
});

smoke('QNumInput: typing a qualified number packs it into the value signal', async () => {
  const input = mount(new QNumInput({label: 'IC50'}));
  assert.equal(input.value.value, null);
  assert.equal(editor(input).value, '');

  type(input, '>10');
  assert.equal(QNum.getQ(input.value.value), QNum.GREATER);
  assert.equal(QNum.getValue(input.value.value), 10);
  assert.equal(input.validity.value, null);

  type(input, '<5.2');
  assert.equal(QNum.getQ(input.value.value), QNum.LESS);
  assert.equal(QNum.getValue(input.value.value), 5.199999999999999);

  type(input, '');
  assert.equal(input.value.value, null, 'blank text is null, as in Dart (qnum_input.dart:32)');
  assert.equal(input.validity.value, null);
  input.dispose();
});

smoke('QNumInput: unparseable text stays on screen and marks the input invalid', async () => {
  const input = mount(new QNumInput({label: 'IC50', value: QNum.create(1, QNum.LESS)}));
  assert.equal(editor(input).value, '<1');

  type(input, '<abc');
  assert.equal(input.validity.value, 'Not a qualified number');
  assert.equal(editor(input).value, '<abc', 'the text is not taken away');
  assert.equal(QNum.getValue(input.value.value), 1, 'the value keeps the last good number');

  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '<abc', 'blur does not overwrite text that is still invalid');

  type(input, '<2');
  assert.equal(input.validity.value, null);
  assert.equal(QNum.getValue(input.value.value), 2);
  input.dispose();
});

smoke('QNumInput: blur redraws the committed value, a write from code redraws it too', async () => {
  const changes = [];
  const input = mount(new QNumInput({label: 'IC50', onChanged: (v) => changes.push(QNum.format(v))}));
  type(input, '> 10.');
  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '>10');

  input.value.value = QNum.create(0.25, QNum.LESS);
  assert.equal(editor(input).value, '<0.25');
  assert.deepEqual(changes, ['>10', '<0.25']);
  input.dispose();
});

smoke('QNumInput: format drives the display, the value stays the packed double', async () => {
  const input = mount(new QNumInput({label: 'IC50', format: (v) => v.toFixed(1),
    value: QNum.create(5.2, QNum.LESS)}));
  assert.equal(editor(input).value, '<5.2');
  assert.equal(QNum.getQ(input.value.value), QNum.LESS);

  type(input, '>1.25');
  assert.equal(QNum.getValue(input.value.value), 1.25, 'the typed number is not rounded');
  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '>1.3', 'only the display is');
  input.dispose();
});

smoke('QNumInput: two editors on one signal, and dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const first = mount(new QNumInput({label: 'From'}));
  const second = mount(new QNumInput({label: 'Also from', bind: first.value}));

  type(first, '<3');
  assert.equal(editor(second).value, '<3');

  second.dispose();
  first.dispose();
  assert.equal(Scope.liveCount, before);

  type(first, '>9');
  assert.equal(QNum.getValue(first.value.value), 3, 'listeners died with the scope');
});
