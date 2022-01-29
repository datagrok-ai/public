export function tryParseJson(s: string): any {
  try {
    return JSON.parse(s);
  } catch (_) {
    return null;
  }
}

export class StringUtils {
  static getUniqueName(initialName: string, existingNames: string[]): string {
    if (!existingNames.includes(initialName))
      return initialName;
    else {
      let counter: number = 1;
      let newName: string = (' ' + initialName + '_' + counter).slice(1);
      while (existingNames.includes(newName)) {
        counter++;
        newName = (' ' + initialName + '_' + counter).slice(1);
      }

      return newName;
    }
  }
}