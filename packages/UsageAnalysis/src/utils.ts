export const colors = {'passed': '#3CB173', 'failed': '#EB6767', 'skipped': '#FFA24A'};

export function getTime(date: Date, format: string = 'en-GB'): string {
  return date.toLocaleString(format, {hour12: false, timeZone: 'GMT'}).replace(',', '');
}

export function getDate(date: Date): string {
  return date.toLocaleDateString('en-US', {month: '2-digit', day: '2-digit', year: 'numeric'});
}
