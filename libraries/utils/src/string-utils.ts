const api = <any>window;

export function tryParseJson(s: string): any {
  try {
    return JSON.parse(s);
  } catch (_) {
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

  static isEmpty(s: string | null): boolean {
    return s == null || s == '';
  }

  static formatNumber(x: number) {
    return _format.format(x);
  }
}

