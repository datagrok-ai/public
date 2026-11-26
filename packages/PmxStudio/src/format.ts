// (a leaf: the dosing schedules format their totals with it, and import nothing)

/** A measured number as it is read: three significant figures (or `digits`), thousands separated ("12,300"), with
 *  no float noise (a float column's 0.3 reads back as 0.30000001192092896: "0.3"); not a number: "—". Every exposure,
 *  AUC, Cmax, concentration, dose and fraction the app prints goes through it. */
export function sig(x: number, digits = 3): string {
  if (typeof x !== 'number' || !Number.isFinite(x))
    return '—';
  const r = Number(x.toPrecision(digits));
  return Math.abs(r) >= 1000 ? r.toLocaleString('en-US', {maximumFractionDigits: 20}) : String(r);
}

/** A value as entered (a dose): every digit kept (3.125, 112.5), only a 32-bit column's float noise off, thousands
 *  separated. sig rounds what was measured or derived; never what the user typed. */
export const entered = (x: number) => sig(x, 7);
