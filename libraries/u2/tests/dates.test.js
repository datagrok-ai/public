/* Dates: the one date-only predicate and the UTC parts everything that shows or edits a UTC date
   goes through. The assertions are written so they are TZ-independent but only BITE away from
   UTC — the suite is run under `TZ=America/New_York`, which is the zone where a UTC midnight read
   locally is the day before. */

// the zone the titles name, set here so the suite bites wherever it is run
process.env.TZ = 'America/New_York';

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {Scope} from '../src/index.js';
import {Dates} from '../src/core/dates.js';

function smoke(name, body) {
  test(name, () => {
    const live = Scope.liveCount;
    body();
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

smoke('isDateOnly: midnight UTC is a date, a local midnight away from UTC is not (TZ=America/New_York)', () => {
  assert.equal(Dates.isDateOnly(new Date('2026-10-02T00:00:00Z')), true);
  assert.equal(Dates.isDateOnly('2026-10-02T00:00:00Z'), true);
  assert.equal(Dates.isDateOnly(Date.UTC(2026, 9, 2)), true);

  // 00:00 LOCAL: a date only where the local zone IS UTC
  const localMidnight = new Date(2026, 9, 2);
  assert.equal(Dates.isDateOnly(localMidnight), localMidnight.getTimezoneOffset() === 0);

  // the old `timestamp()` predicate tested UTC h/m/s and ignored the milliseconds
  assert.equal(Dates.isDateOnly(new Date('2026-10-02T00:00:00.500Z')), false);
  assert.equal(Dates.isDateOnly(new Date('2026-10-02T12:00:00Z')), false);
  assert.equal(Dates.isDateOnly(new Date('nonsense')), false);
});

smoke('isDateOnlyMicros: the null sentinel and a non-finite value say nothing either way', () => {
  const FLOAT_NULL = 2.6789344063684636e-34;
  assert.equal(Dates.isDateOnlyMicros(FLOAT_NULL), true);
  assert.equal(Dates.isDateOnlyMicros(NaN), true);
  assert.equal(Dates.isDateOnlyMicros(Infinity), true);

  const day = Date.UTC(2026, 9, 2) * 1000;
  assert.equal(Dates.isDateOnlyMicros(day), true);
  assert.equal(Dates.isDateOnlyMicros(day + 12 * 3600 * 1e6), false);
  assert.equal(Dates.isDateOnlyMicros(0), true);
  assert.equal(Dates.MICROSECONDS_PER_DAY, 86400000000);
});

smoke('utcParts / fromUtcParts round-trip across a DST boundary (TZ=America/New_York)', () => {
  // the US falls back on 2026-11-01 and springs forward on 2026-03-08
  for (const iso of ['2026-03-07T23:30:00Z', '2026-03-08T07:30:00Z', '2026-11-01T05:30:00Z',
    '2026-11-01T06:30:00Z', '2026-06-15T00:00:00Z']) {
    const date = new Date(iso);
    const parts = Dates.utcParts(date);
    assert.equal(Dates.fromUtcParts(parts).toISOString(), iso.replace('Z', '.000Z'), iso);
  }

  const parts = Dates.utcParts(new Date('2026-01-02T03:04:00Z'));
  assert.deepEqual(parts, {year: 2026, month: 1, day: 2, hours: 3, minutes: 4});
  assert.equal(Dates.fromUtcParts({year: 2026, month: 1, day: 2}).toISOString(), '2026-01-02T00:00:00.000Z');
});
