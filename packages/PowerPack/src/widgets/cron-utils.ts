export interface CronFields {
  minute: string;
  hour: string;
  dayOfMonth: string;
  month: string;
  dayOfWeek: string;
}

export interface CronPreset {
  label: string;
  value: string;
}

export interface CronFieldMeta {
  key: keyof CronFields;
  label: string;
  min: number;
  max: number;
}

export const CRON_FIELD_META: CronFieldMeta[] = [
  {key: 'minute', label: 'Minute', min: 0, max: 59},
  {key: 'hour', label: 'Hour', min: 0, max: 23},
  {key: 'dayOfMonth', label: 'Day of month', min: 1, max: 31},
  {key: 'month', label: 'Month', min: 1, max: 12},
  {key: 'dayOfWeek', label: 'Day of week', min: 0, max: 7},
];

export const CRON_PRESETS: CronPreset[] = [
  {label: 'Every minute', value: '* * * * *'},
  {label: 'Every 5 minutes', value: '*/5 * * * *'},
  {label: 'Every 15 minutes', value: '*/15 * * * *'},
  {label: 'Every 30 minutes', value: '*/30 * * * *'},
  {label: 'Hourly', value: '0 * * * *'},
  {label: 'Daily at midnight', value: '0 0 * * *'},
  {label: 'Daily at 9:00 AM', value: '0 9 * * *'},
  {label: 'Weekly on Monday', value: '0 0 * * 1'},
  {label: 'Monthly on the 1st', value: '0 0 1 * *'},
];

export function parseCron(expr: string): CronFields | null {
  if (!expr)
    return null;
  const parts = expr.trim().split(/\s+/);
  if (parts.length !== 5)
    return null;
  return {
    minute: parts[0],
    hour: parts[1],
    dayOfMonth: parts[2],
    month: parts[3],
    dayOfWeek: parts[4],
  };
}

export function cronToString(fields: CronFields): string {
  return `${fields.minute} ${fields.hour} ${fields.dayOfMonth} ${fields.month} ${fields.dayOfWeek}`;
}

function isValidToken(token: string, min: number, max: number): boolean {
  if (token === '*')
    return true;

  for (const part of token.split(',')) {
    const stepParts = part.split('/');
    if (stepParts.length > 2)
      return false;

    const rangePart = stepParts[0];
    const step = stepParts[1];

    if (step !== undefined) {
      const s = parseInt(step);
      if (isNaN(s) || s < 1)
        return false;
    }

    if (rangePart === '*')
      continue;

    const rangeBounds = rangePart.split('-');
    if (rangeBounds.length > 2)
      return false;

    for (const b of rangeBounds) {
      const n = parseInt(b);
      if (isNaN(n) || n < min || n > max)
        return false;
    }

    if (rangeBounds.length === 2) {
      const lo = parseInt(rangeBounds[0]);
      const hi = parseInt(rangeBounds[1]);
      if (lo > hi)
        return false;
    }
  }
  return true;
}

export function isValidCron(expr: string): boolean {
  const fields = parseCron(expr);
  if (!fields)
    return false;
  for (const meta of CRON_FIELD_META) {
    if (!isValidToken(fields[meta.key], meta.min, meta.max))
      return false;
  }
  return true;
}

/** Returns first validation error, or null if valid. */
export function validateCron(expr: string): string | null {
  const fields = parseCron(expr);
  if (!fields)
    return 'Expression must have exactly 5 fields';
  for (const meta of CRON_FIELD_META) {
    if (!isValidToken(fields[meta.key], meta.min, meta.max))
      return `${meta.label} must be between ${meta.min}\u2013${meta.max}`;
  }
  return null;
}

/** Returns per-field error messages for invalid fields. */
export function validateCronPerField(fields: CronFields): {[key: string]: string} {
  const errors: {[key: string]: string} = {};
  for (const meta of CRON_FIELD_META) {
    if (!isValidToken(fields[meta.key], meta.min, meta.max))
      errors[meta.key] = `${meta.label} must be between ${meta.min}\u2013${meta.max}`;
  }
  return errors;
}

export function findMatchingPreset(expr: string): CronPreset | null {
  const normalized = expr.trim().replace(/\s+/g, ' ');
  for (const preset of CRON_PRESETS) {
    if (preset.value === normalized)
      return preset;
  }
  return null;
}

// --- Next run computation ---

function expandField(expr: string, min: number, max: number): Set<number> {
  const values = new Set<number>();
  for (const part of expr.split(',')) {
    const stepParts = part.split('/');
    const rangePart = stepParts[0];
    const step = stepParts[1] ? parseInt(stepParts[1]) : 1;

    let lo: number, hi: number;
    if (rangePart === '*') {
      lo = min;
      hi = max;
    }
    else if (rangePart.includes('-')) {
      const bounds = rangePart.split('-');
      lo = parseInt(bounds[0]);
      hi = parseInt(bounds[1]);
    }
    else {
      lo = parseInt(rangePart);
      hi = stepParts[1] ? max : lo;
    }

    for (let v = lo; v <= hi; v += step)
      values.add(v);
  }
  return values;
}

export function getNextRun(expr: string, from?: Date): Date | null {
  const fields = parseCron(expr);
  if (!fields || !isValidCron(expr))
    return null;

  const minuteSet = expandField(fields.minute, 0, 59);
  const hourSet = expandField(fields.hour, 0, 23);
  const domSet = expandField(fields.dayOfMonth, 1, 31);
  const monthSet = expandField(fields.month, 1, 12);
  const dowSet = expandField(fields.dayOfWeek, 0, 7);
  if (dowSet.has(7)) dowSet.add(0);
  if (dowSet.has(0)) dowSet.add(7);

  const d = new Date(from || new Date());
  d.setSeconds(0, 0);
  d.setMinutes(d.getMinutes() + 1);

  const limit = new Date(d);
  limit.setFullYear(limit.getFullYear() + 2);

  while (d < limit) {
    if (!monthSet.has(d.getMonth() + 1)) {
      d.setMonth(d.getMonth() + 1, 1);
      d.setHours(0, 0, 0, 0);
      continue;
    }
    if (!domSet.has(d.getDate()) || !dowSet.has(d.getDay())) {
      d.setDate(d.getDate() + 1);
      d.setHours(0, 0, 0, 0);
      continue;
    }
    if (!hourSet.has(d.getHours())) {
      d.setHours(d.getHours() + 1, 0, 0, 0);
      continue;
    }
    if (!minuteSet.has(d.getMinutes())) {
      d.setMinutes(d.getMinutes() + 1);
      continue;
    }
    return new Date(d);
  }
  return null;
}

// --- Human-readable description ---

const DAY_NAMES = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'];
const DAY_SHORT = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'];
const MONTH_NAMES = ['', 'January', 'February', 'March', 'April', 'May', 'June',
  'July', 'August', 'September', 'October', 'November', 'December'];

function fmtTime(hour: string, minute: string): string {
  const h = parseInt(hour);
  const m = parseInt(minute);
  return `${h < 10 ? '0' + h : h}:${m < 10 ? '0' + m : m}`;
}

function ordinal(n: number): string {
  return n === 1 ? '1st' : n === 2 ? '2nd' : n === 3 ? '3rd' : n + 'th';
}

export function cronToDescription(expr: string): string {
  const fields = parseCron(expr);
  if (!fields)
    return 'Invalid expression';

  const {minute, hour, dayOfMonth, month, dayOfWeek} = fields;

  if (minute === '*' && hour === '*' && dayOfMonth === '*' && month === '*' && dayOfWeek === '*')
    return 'Every minute';

  if (minute.startsWith('*/') && hour === '*' && dayOfMonth === '*' && month === '*' && dayOfWeek === '*') {
    const step = minute.slice(2);
    return step === '1' ? 'Every minute' : `Every ${step} minutes`;
  }

  if (!minute.includes('*') && !minute.includes('/') && hour === '*' &&
      dayOfMonth === '*' && month === '*' && dayOfWeek === '*') {
    const m = parseInt(minute);
    return m === 0 ? 'Every hour' : `At minute ${m} of every hour`;
  }

  if (minute === '0' && hour.startsWith('*/') && dayOfMonth === '*' && month === '*' && dayOfWeek === '*') {
    const step = hour.slice(2);
    return step === '1' ? 'Every hour' : `Every ${step} hours`;
  }

  if (!minute.includes('*') && !minute.includes('/') &&
      !hour.includes('*') && !hour.includes('/') &&
      dayOfMonth === '*' && month === '*' && dayOfWeek === '*')
    return `Daily at ${fmtTime(hour, minute)}`;

  if (!minute.includes('*') && !minute.includes('/') &&
      !hour.includes('*') && !hour.includes('/') &&
      dayOfMonth === '*' && month === '*' &&
      !dayOfWeek.includes('*') && !dayOfWeek.includes('/')) {
    const days = dayOfWeek.split(',').map((d) => {
      const n = parseInt(d);
      return (n >= 0 && n <= 7) ? DAY_NAMES[n] : d;
    });
    return `At ${fmtTime(hour, minute)} on ${days.join(', ')}`;
  }

  if (!minute.includes('*') && !minute.includes('/') &&
      !hour.includes('*') && !hour.includes('/') &&
      !dayOfMonth.includes('*') && !dayOfMonth.includes('/') &&
      month === '*' && dayOfWeek === '*')
    return `At ${fmtTime(hour, minute)} on the ${ordinal(parseInt(dayOfMonth))} of every month`;

  if (!minute.includes('*') && !minute.includes('/') &&
      !hour.includes('*') && !hour.includes('/') &&
      !dayOfMonth.includes('*') && !dayOfMonth.includes('/') &&
      !month.includes('*') && !month.includes('/') && dayOfWeek === '*') {
    const mon = parseInt(month);
    const monthName = (mon >= 1 && mon <= 12) ? MONTH_NAMES[mon] : month;
    return `At ${fmtTime(hour, minute)} on ${monthName} ${ordinal(parseInt(dayOfMonth))}`;
  }

  return expr;
}

/** Translate a weekday token to human-readable short names. Returns null for wildcards. */
export function dowTokenToNames(token: string): string | null {
  if (token === '*')
    return null;

  const parts = token.split(',');
  const names: string[] = [];
  for (const part of parts) {
    if (part.includes('/'))
      return null;
    if (part.includes('-')) {
      const bounds = part.split('-');
      const lo = parseInt(bounds[0]);
      const hi = parseInt(bounds[1]);
      if (isNaN(lo) || isNaN(hi) || lo < 0 || lo > 7 || hi < 0 || hi > 7)
        return null;
      names.push(`${DAY_SHORT[lo]}\u2013${DAY_SHORT[hi]}`);
    }
    else {
      const n = parseInt(part);
      if (isNaN(n) || n < 0 || n > 7)
        return null;
      names.push(DAY_SHORT[n]);
    }
  }
  return names.join(', ');
}

const SHORT_MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

export function formatNextRun(date: Date): string {
  const month = SHORT_MONTHS[date.getMonth()];
  const day = date.getDate();
  const year = date.getFullYear();
  const h = date.getHours();
  const m = date.getMinutes();
  const hh = h < 10 ? '0' + h : '' + h;
  const mm = m < 10 ? '0' + m : '' + m;
  return `${month} ${day}, ${year} at ${hh}:${mm}`;
}
