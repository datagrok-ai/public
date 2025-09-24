/* eslint-disable max-len */
type Nullable<T> = T | null;


export abstract class Matcher {
  abstract match(x: any): boolean;
  abstract toSql(variable: string): string;
}

export class StringInListMatcher extends Matcher {
  values: string[];

  constructor(values: string[]) {
    super();
    this.values = values;
  }

  match(x: string | null): boolean {
    if (x === null || x === undefined)
      return false;
    return this.values.includes(x);
  }

  toSql(variable: string): string {
    if (this.values.length === 0)
      return '(1 = 1)';

    const escapedValues = this.values
      .map((v) => `'${v.replace(/'/g, '\'\'')}'`)
      .join(', ');
    return `(${variable} IN (${escapedValues}))`;
  }
}

export class NumericMatcher {
  static readonly NONE = 'none';
  static readonly EQUALS = '=';
  static readonly NOT_EQUALS = '!=';
  static readonly GT = '>';
  static readonly GTOE = '>=';
  static readonly LT = '<';
  static readonly LTOE = '<=';
  static readonly IN = 'in';
  static readonly NOT_IN = 'not in';
  static readonly RANGE = '-';
  static readonly IS_NULL = 'is null';
  static readonly IS_NOT_NULL = 'is not null';

  static readonly _num = '[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?';
  static readonly numRegex = new RegExp(`^${NumericMatcher._num}`);
  static readonly unaryRegex = /(?:[$]{.*})?\s*(!=|>=|<=|=|>|<|\bnot in\b|\bin\b)+/;
  static readonly rangeRegex = new RegExp(`(${NumericMatcher._num})\\s?(?:-|\\.\\.)\\s?(${NumericMatcher._num})`);
  static readonly inRegex = /\((.*)\)/;
  static readonly isNullRegExp = /(^is null|^is not null|\sis null|\sis not null)$/;

  static readonly unaryOperators = [
    NumericMatcher.EQUALS, NumericMatcher.NOT_EQUALS,
    NumericMatcher.GT, NumericMatcher.GTOE,
    NumericMatcher.LT, NumericMatcher.LTOE,
    NumericMatcher.IN, NumericMatcher.NOT_IN
  ];

  static recentMatchers = new Map<string, NumericMatcher>();

  op?: string;
  v1?: number;
  v2?: number;
  values: number[] = [];
  expression?: string;
  allowedPatterns?: string[];

  constructor(op?: string, v1?: number, v2?: number) {
    this.op = op;
    this.v1 = v1;
    this.v2 = v2;
    if (v1 !== undefined) this.values.push(v1);
    if (v2 !== undefined) this.values.push(v2);
  }

  static parse(query: string, options?: {
    orElse?: () => NumericMatcher,
    matchEmptyPattern?: boolean
  }): Nullable<NumericMatcher> {
    const {orElse, matchEmptyPattern} = options || {};
    const fail = () => orElse ? orElse() : null;

    if (matchEmptyPattern && (!query || query.trim() === ''))
      return new NumericMatcher(NumericMatcher.NONE);

    const x = parseFloat(query);
    if (!isNaN(x))
      return new NumericMatcher(NumericMatcher.EQUALS, x);

    const match = NumericMatcher.unaryRegex.exec(query);
    if (match && NumericMatcher.unaryOperators.includes(match[1])) {
      const op = match[1];
      const matcher = new NumericMatcher(op);
      matcher.expression = query;
      const rest = query.substring(match.index + match[0].length).trim();

      if (op === NumericMatcher.IN || op === NumericMatcher.NOT_IN) {
        const matchValues = NumericMatcher.inRegex.exec(rest);
        if (matchValues) {
          const values = matchValues[1].split(',').map((s) => parseFloat(s.trim()));
          if (values.some(isNaN)) return fail();
          matcher.values = values;
          return matcher;
        }
      } else {
        const matchValue1 = NumericMatcher.numRegex.exec(rest.toLowerCase());
        if (matchValue1) {
          matcher.v1 = parseFloat(matchValue1[0]);
          matcher.values = [matcher.v1];
          return matcher;
        }
      }
    }

    const upQuery = query.toUpperCase();
    if (['NONE', 'EMPTY', 'NA', 'N/A'].includes(upQuery))
      return new NumericMatcher(NumericMatcher.NONE);

    const rangeMatch = NumericMatcher.rangeRegex.exec(query);
    if (rangeMatch) {
      return new NumericMatcher(
        NumericMatcher.RANGE,
        parseFloat(rangeMatch[1]),
        parseFloat(rangeMatch[2])
      );
    }

    const nullMatch = NumericMatcher.isNullRegExp.exec(query);
    if (nullMatch)
      return new NumericMatcher(nullMatch[1].trim());

    return fail();
  }

  static matches(query: string, x: number, matchEmptyPattern = false): boolean {
    let matcher: NumericMatcher | undefined | null = this.recentMatchers.get(query);
    if (!matcher) {
      matcher = this.parse(query, {matchEmptyPattern});
      if (matcher)
        this.recentMatchers.set(query, matcher);
    }
    return matcher?.match(x) === true;
  }

  match(x?: number | null): boolean {
    if (x === null || x === undefined)
      return this.op === NumericMatcher.NONE || this.op === NumericMatcher.IS_NULL;

    switch ((this.op || '').trim()) {
      case NumericMatcher.NONE: return false; // we already know x != null;
      case NumericMatcher.GT: return x > this.v1!;
      case NumericMatcher.GTOE: return x >= this.v1!;
      case NumericMatcher.LT: return x < this.v1!;
      case NumericMatcher.LTOE: return x <= this.v1!;
      case NumericMatcher.EQUALS: return x === this.v1!;
      case NumericMatcher.NOT_EQUALS: return x !== this.v1!;
      case NumericMatcher.IN: return this.values.includes(x);
      case NumericMatcher.NOT_IN: return !this.values.includes(x);
      case NumericMatcher.RANGE: return x >= this.v1! && x <= this.v2!;
      case NumericMatcher.IS_NULL: return false;
      case NumericMatcher.IS_NOT_NULL: return true;
      default: throw new Error(`Unknown operation "${this.op}"`);
    }
  }

  validatePattern(): string | null {
    if ((!this.allowedPatterns || this.allowedPatterns.includes(this.op!)) &&
      (NumericMatcher.unaryOperators.includes(this.op!) ||
        [NumericMatcher.NONE, NumericMatcher.RANGE, NumericMatcher.IS_NULL, NumericMatcher.IS_NOT_NULL].includes(this.op!)))
      return null;
    return `Unknown operation "${this.op}"`;
  }

  toSql(variable: string): string {
    switch (this.op) {
      case NumericMatcher.RANGE:
        return `(${variable} >= ${this.v1} AND ${variable} <= ${this.v2})`;
      case NumericMatcher.IN:
        return `(${variable} IN (${this.values.join(',')}))`;
      case NumericMatcher.NONE:
        return `(1 = 1)`;
      case NumericMatcher.IS_NULL:
      case NumericMatcher.IS_NOT_NULL:
        return `(${variable} ${this.op})`;
      default:
        return `(${variable} ${this.op} ${this.v1})`;
    }
  }
}

export class StringMatcher extends Matcher {
  value: string;

  constructor(value: string) {
    super();
    this.value = value;
  }

  match(x: string | null): boolean {
    if (x === null || x === undefined)
      return false;
    return x.toLowerCase().includes(this.value.toLowerCase());
  }

  toSql(variable: string): string {
    // Escape single quotes for SQL, then prepare for a case-insensitive LIKE search
    const escapedValue = this.value.replace(/'/g, '\'\'');
    return `(${variable} ILIKE '%${escapedValue}%')`;
  }
}
