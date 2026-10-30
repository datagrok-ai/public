const api: any = (typeof window !== 'undefined' ? window : global.window) as any;

export function tryParseJson(s: string, options?: {warning?: string, reviver?: (key: string, value: any) => any}): any {
  try {
    return JSON.parse(s, options?.reviver);
  } catch (e) {
    if (options?.warning)
      console.warn(`${options.warning}: malformed JSON: ${s}`, e);
    return null;
  }
}

const _format = new Intl.NumberFormat('en-US', {minimumFractionDigits: 2, maximumFractionDigits: 2});

export class StringUtils {


  /**
   * Returns a name that does not exist in [existingNames].
   * If [existingNames] does not contain [initialName], returns [initialName].
   * Otherwise, tries [choices], and if the names are taken already, returns a string in a form of 'initialName (i)'.
   */
  static getUniqueName(initialName: string, existingNames: string[],
    options?: { auto?: boolean, idx?: number, render?: Function, choices?: string[] }): string {
    return api.grok_Utils_GetUniqueName(
      initialName, existingNames, options?.auto, options?.idx, options?.render, options?.choices);
  }

  /** [s] as a regular expression that matches the text itself. */
  static escapeRegExp(s: string): string {
    return s.replace(/[.*+?^${}()|[\]\\\/-]/g, '\\$&');
  }

  static isEmpty(s: string | null): boolean {
    return s == null || s == '';
  }

  static formatNumber(x: number) {
    return _format.format(x);
  }
}

