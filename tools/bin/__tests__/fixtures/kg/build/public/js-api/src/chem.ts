export namespace chem {
  export interface Options {
    timeout: number;
  }

  export interface SearchOptions extends Options {
    limit: number;
  }

  export function similarity(a: string, b: string): number {
    return 0;
  }
}
