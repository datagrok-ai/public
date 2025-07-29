import {ACTION_KEYS, PRINTABLE_LEN} from './constants';

export function debounce<T extends(...args: any[]) => void>(
  func: T,
  delay: number,
): (...args: Parameters<T>) => void {
  let timeoutId: ReturnType<typeof setTimeout>;

  return function(...args: Parameters<T>) {
    clearTimeout(timeoutId);
    timeoutId = setTimeout(() => {
      func(...args);
    }, delay);
  };
}

export function isActionKey(key: string): boolean {
  return (key.length === PRINTABLE_LEN) || ACTION_KEYS.has(key);
}
