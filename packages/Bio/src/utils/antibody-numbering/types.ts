/** Position → amino-acid entry in the numbering detail JSON. */
export interface ImmunumNumberingEntry {
  /** Position code (e.g. "1", "27A"). */
  position: string;
  /** Single-letter residue at this position. */
  aa: string;
}

/** Per-row numbering result produced by the immunum worker. */
export interface ImmunumNumberingRow {
  /** Comma-separated position codes in scheme order. Empty string on failure. */
  positionNames: string;
  /** 'Heavy' / 'Light' / '' — UI-facing chain group. */
  chainType: string;
  /** Immunum chain code: H, K, L, A, B, G, D, or ''. */
  chainCode: string;
  /** Ordered list of {position, aa} entries for the aligned region. */
  numberingDetail: ImmunumNumberingEntry[];
  /** Map from position code → character index in the extracted (gap-free) sequence. */
  numberingMap: Record<string, number>;
  /** Alignment confidence in [0, 1]. */
  confidence: number;
  /** Error message, empty string on success. */
  error: string;
}

export type ImmunumWorkerRequest =
  | {op: 'init'}
  | {
      op: 'number';
      sequences: string[];
      /** Immunum scheme: 'imgt' or 'kabat'. */
      scheme: string;
      /** Case-insensitive chain codes (e.g. ['H', 'K', 'L']). Defaults to H/K/L. */
      chains?: string[];
      /** Minimum confidence in [0, 1]; null/undefined → library default (0.5). */
      minConfidence?: number | null;
    };

export type ImmunumWorkerResponse =
  | {ok: true; rows?: ImmunumNumberingRow[]}
  | {ok: false; error: string};

/** Supported schemes by immunum. Chothia/AHo are not supported — the UI engine
 *  dropdown lists all schemes globally so the immunum path falls back for the
 *  unsupported ones by returning empty rows with errors. */
export const IMMUNUM_SCHEMES = ['imgt', 'kabat'] as const;
export type ImmunumScheme = typeof IMMUNUM_SCHEMES[number];
