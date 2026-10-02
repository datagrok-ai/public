/** A table: columns, rows, selection and filter. */
export class DataFrame {
  private _name: string = '';

  /** Name of the table. */
  get name(): string {
    return this._name;
  }

  set name(s: string) {
    this._name = s;
  }

  get rowCount(): number {
    return 0;
  }

  /** @deprecated Use fromCsv. */
  static fromText(csv: string): DataFrame {
    return new DataFrame();
  }

  static fromCsv(csv: string, options?: {delimiter?: string}): DataFrame {
    return new DataFrame();
  }

  protected hidden(): void {
  }
}

export interface IDisposable {
  dispose(): void;
}

export enum LogLevel {
  Info = 'info',
  Error = 'error',
}

export type Predicate = (x: number) => boolean;

class Internal {
  m(): void {
  }
}
