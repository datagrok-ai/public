/**
 * AnnotationRenderer — overlays annotation visuals on top of sequence cells.
 *
 * Performance contract:
 * - When no annotations exist, `hasAnnotations()` returns false and callers skip all drawing.
 * - Parsed annotation data is cached by column version; re-parsed only on change.
 * - The companion row-data column (hidden, name starting with ~) is re-looked up when the
 *   sequence column version changes (so it's found even when created after first render).
 */
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from './macromolecule/consts';
import {
  SeqAnnotation, SeqAnnotationHit, RowAnnotationData,
  AnnotationCategory, REGION_BG_OPACITY, LIABILITY_UNDERLINE_HEIGHT,
  ANNOTATION_COLORS,
} from './macromolecule/annotations';

/** Cached annotation state for a column. */
interface AnnotationCache {
  colVersion: number;
  annotations: SeqAnnotation[];
  /** Map from 0-based position index → region annotation color */
  positionRegionColors: Map<number, string>;
  /** region annotation names by position (for tooltip) */
  positionRegionNames: Map<number, string>;
}

interface RowDataCache {
  colVersion: number;
  data: (RowAnnotationData | null)[];
}

export class AnnotationRenderer {
  private _cache: AnnotationCache | null = null;
  private _rowCache: RowDataCache | null = null;
  private _annotCol: DG.Column<string> | null = null;
  /** Version of tableCol at which _annotCol was last looked up.
   *  Re-lookup when version changes (setting .annotationColumnName tag bumps it). */
  private _annotColLookupVersion: number = -1;

  constructor(
    private tableCol: DG.Column<string>,
  ) {}

  /** Fast check — returns false when no annotations exist, so the render loop can skip entirely. */
  hasAnnotations(): boolean {
    this._ensureCache();
    return this._cache !== null && this._cache.annotations.length > 0;
  }

  /** Returns region annotations (column-level). */
  getAnnotations(): SeqAnnotation[] {
    this._ensureCache();
    return this._cache?.annotations ?? [];
  }

  /** Returns the region name at a given position for a specific row.
   *  Checks per-row region spans first; falls back to column-level only when
   *  the row has no per-row region data at all. */
  getRegionNameAtPosition(posIdx: number, rowIdx?: number): string | null {
    // Check per-row region spans first
    if (rowIdx != null) {
      const hit = this._findRowRegionHit(rowIdx, posIdx);
      if (hit) {
        const annot = this._cache?.annotations.find((a) => a.id === hit.annotationId);
        return annot?.name ?? null;
      }
      // If this row has per-row region data, don't fall back to column-level
      if (this._rowHasRegionData(rowIdx))
        return null;
    }
    // Fall back to column-level position map
    this._ensureCache();
    return this._cache?.positionRegionNames.get(posIdx) ?? null;
  }

  /** Returns liability hits (non-region) for a specific row and position index. */
  getHitsAtPosition(rowIdx: number, posIdx: number): SeqAnnotationHit[] {
    this._ensureRowCache();
    if (!this._rowCache) return [];
    const rowData = this._rowCache.data[rowIdx];
    if (!rowData) return [];
    return rowData.filter((h) => {
      if (h.endPositionIndex != null) return false; // Skip region spans
      return h.positionIndex === posIdx || (h.matchedMonomers.length > 1 && posIdx > h.positionIndex && posIdx < h.positionIndex + h.matchedMonomers.length);
    });
  }

  /** Draws annotation background + underline for a single position in a cell.
   *  Called from MonomerPlacer.render() for each visible position.
   *  @param textBottomY If provided, underline is drawn just below this y-coordinate
   *    (use in single-line mode so the underline sticks to the letters). */
  drawPositionBackground(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    posIdx: number, rowIdx: number, textBottomY?: number,
  ): void {
    this._ensureCache();
    if (!this._cache) return;

    this._ensureRowCache();

    // Region background — check per-row region spans first, then column-level.
    // If per-row region data exists for this row, use ONLY that (don't fall back
    // to column-level, which has generic positions that don't match this row).
    let regionColor: string | undefined;
    const rowRegionHit = this._findRowRegionHit(rowIdx, posIdx);
    if (rowRegionHit) {
      const annot = this._cache.annotations.find((a) => a.id === rowRegionHit.annotationId);
      regionColor = annot?.color ?? (annot ? this._getDefaultRegionColor(annot) : undefined);
    }
    if (!regionColor && !this._rowHasRegionData(rowIdx))
      regionColor = this._cache.positionRegionColors.get(posIdx);

    if (regionColor) {
      g.save();
      g.globalAlpha = REGION_BG_OPACITY;
      g.fillStyle = regionColor;
      g.fillRect(x, y, w, h);
      g.restore();
    }

    // Liability underline — snapped to text bottom when textBottomY is provided
    if (this._rowCache) {
      const rowData = this._rowCache.data[rowIdx];
      if (rowData) {
        for (const hit of rowData) {
          if (hit.endPositionIndex != null) continue; // Skip region spans
          if (posIdx >= hit.positionIndex && posIdx < hit.positionIndex + hit.matchedMonomers.length) {
            const annot = this._cache.annotations.find((a) => a.id === hit.annotationId);
            if (annot?.color) {
              const underlineY = textBottomY != null
                ? textBottomY + 1
                : y + h - LIABILITY_UNDERLINE_HEIGHT;
              g.fillStyle = annot.color;
              g.fillRect(x, underlineY, w, LIABILITY_UNDERLINE_HEIGHT);
            }
            break; // Only draw one underline
          }
        }
      }
    }
  }

  /** Builds a tooltip HTML snippet for annotations at a given position+row. */
  getTooltipInfo(posIdx: number, rowIdx: number): string[] {
    this._ensureCache();
    const lines: string[] = [];
    if (!this._cache) return lines;

    this._ensureRowCache();

    // Region info — check per-row region spans first; fall back to column-level
    // only when the row has no per-row region data at all.
    let regionName: string | null = null;
    const rowRegionHit = this._findRowRegionHit(rowIdx, posIdx);
    if (rowRegionHit) {
      const annot = this._cache.annotations.find((a) => a.id === rowRegionHit.annotationId);
      if (annot) {
        regionName = annot.name;
        const scheme = annot.sourceScheme ?? '';
        lines.push(`Region: ${annot.name}${scheme ? ` (${scheme} ${annot.start}-${annot.end})` : ''}`);
      }
    }
    if (!regionName && !this._rowHasRegionData(rowIdx)) {
      const colRegionName = this._cache.positionRegionNames.get(posIdx);
      if (colRegionName) {
        const annot = this._cache.annotations.find((a) => a.name === colRegionName && a.category === AnnotationCategory.Structure);
        const scheme = annot?.sourceScheme ?? '';
        lines.push(`Region: ${colRegionName}${scheme ? ` (${scheme} ${annot?.start}-${annot?.end})` : ''}`);
      }
    }

    // Liability hits
    if (this._rowCache) {
      const rowData = this._rowCache.data[rowIdx];
      if (rowData) {
        for (const hit of rowData) {
          if (hit.endPositionIndex != null) continue; // Skip region spans
          if (posIdx >= hit.positionIndex && posIdx < hit.positionIndex + hit.matchedMonomers.length) {
            const annot = this._cache.annotations.find((a) => a.id === hit.annotationId);
            if (annot) {
              const severityLabel = annot.severity ? ` - ${annot.severity.charAt(0).toUpperCase() + annot.severity.slice(1)}` : '';
              lines.push(`\u26A0 ${annot.name}${severityLabel}`);
            }
          }
        }
      }
    }

    return lines;
  }

  /** Finds a per-row region span hit covering the given position, or null. */
  private _findRowRegionHit(rowIdx: number, posIdx: number): SeqAnnotationHit | null {
    this._ensureRowCache();
    if (!this._rowCache) return null;
    const rowData = this._rowCache.data[rowIdx];
    if (!rowData) return null;
    for (const hit of rowData) {
      if (hit.endPositionIndex != null && posIdx >= hit.positionIndex && posIdx <= hit.endPositionIndex)
        return hit;
    }
    return null;
  }

  /** Returns true if the row has any per-row region span data (hits with endPositionIndex set).
   *  When true, column-level region mapping should NOT be used as fallback for this row. */
  private _rowHasRegionData(rowIdx: number): boolean {
    this._ensureRowCache();
    if (!this._rowCache) return false;
    const rowData = this._rowCache.data[rowIdx];
    if (!rowData) return false;
    return rowData.some((h) => h.endPositionIndex != null);
  }

  // --- Internal caching ---

  private _ensureCache(): void {
    const currentVersion = this.tableCol.version;
    if (this._cache && this._cache.colVersion === currentVersion) return;

    const annotTag = this.tableCol.getTag(bioTAGS.annotations);
    if (!annotTag) {
      this._cache = null;
      return;
    }

    let annotations: SeqAnnotation[];
    try {
      annotations = JSON.parse(annotTag);
    } catch {
      this._cache = null;
      return;
    }

    if (!annotations || annotations.length === 0) {
      this._cache = null;
      return;
    }

    // Build position maps for structure annotations
    const positionRegionColors = new Map<number, string>();
    const positionRegionNames = new Map<number, string>();
    const posList = this._getPosList();

    for (const annot of annotations) {
      if (annot.category !== AnnotationCategory.Structure) continue;
      if (annot.start == null || annot.end == null) continue;

      const startIdx = posList.indexOf(annot.start);
      const endIdx = posList.indexOf(annot.end);
      if (startIdx < 0 || endIdx < 0) continue;

      const color = annot.color ?? this._getDefaultRegionColor(annot);
      for (let i = startIdx; i <= endIdx; i++) {
        positionRegionColors.set(i, color);
        positionRegionNames.set(i, annot.name);
      }
    }

    this._cache = {colVersion: currentVersion, annotations, positionRegionColors, positionRegionNames};
  }

  private _ensureRowCache(): void {
    // Re-lookup companion column when the sequence column version changes.
    // This ensures we pick up the companion column if it was created after
    // the first render (e.g., liability scan creates it later).
    const seqVersion = this.tableCol.version;
    if (this._annotColLookupVersion !== seqVersion) {
      this._annotColLookupVersion = seqVersion;
      const annotColName = this.tableCol.getTag(bioTAGS.annotationColumnName);
      if (annotColName) {
        try {
          this._annotCol = this.tableCol.dataFrame.columns.byName(annotColName) as DG.Column<string> | null;
        } catch {
          this._annotCol = null;
        }
      } else {
        this._annotCol = null;
      }
    }

    if (!this._annotCol) {
      this._rowCache = null;
      return;
    }

    const currentVersion = this._annotCol.version;
    if (this._rowCache && this._rowCache.colVersion === currentVersion) return;

    const data: (RowAnnotationData | null)[] = new Array(this._annotCol.length);
    for (let i = 0; i < this._annotCol.length; i++) {
      const raw = this._annotCol.get(i);
      if (!raw) { data[i] = null; continue; }
      try { data[i] = JSON.parse(raw); } catch { data[i] = null; }
    }

    this._rowCache = {colVersion: currentVersion, data};
  }

  private _getPosList(): string[] {
    const posNamesTag = this.tableCol.getTag(bioTAGS.positionNames);
    if (!posNamesTag) return [];
    return posNamesTag.split(', ');
  }

  private _getDefaultRegionColor(annot: SeqAnnotation): string {
    if (annot.name.startsWith('CDR')) {
      const idx = parseInt(annot.name.replace('CDR', '')) - 1;
      return ANNOTATION_COLORS.structure.CDR[idx % ANNOTATION_COLORS.structure.CDR.length];
    }
    if (annot.name.startsWith('FR')) {
      const idx = parseInt(annot.name.replace('FR', '')) - 1;
      return ANNOTATION_COLORS.structure.FR[idx % ANNOTATION_COLORS.structure.FR.length];
    }
    return '#90CAF9';
  }
}
