export const error = (s: string) => console.log('\x1b[31m%s\x1b[0m', s);
export const info = (s: string) => console.log('\x1b[32m%s\x1b[0m', s);
export const warn = (s: string) => console.log('\x1b[33m%s\x1b[0m', s);

export const success = info;
export const fail = error;
