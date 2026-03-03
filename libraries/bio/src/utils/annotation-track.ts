/**
 * AnnotationTrack — MSA header track displaying colored region bands.
 *
 * Shows FR/CDR regions (or any structure annotations) as colored horizontal
 * rectangles with centered labels. This track is designed to be the primary
 * (first) track in the MSA header when annotation data exists.
 */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {MSAHeaderTrack} from './sequence-position-scroller';
import {TAGS as bioTAGS} from './macromolecule/consts';
import {SeqAnnotation, AnnotationCategory, ANNOTATION_COLORS} from './macromolecule/annotations';

/** Compact track height for the annotation band. */
const ANNOTATION_TRACK_HEIGHT = 20;
const ANNOTATION_TRACK_MIN_HEIGHT = 14;

const ANNOTATION_FONT = 'bold 9px Roboto, Roboto Local, sans-serif';
const ANNOTATION_LABEL_COLOR = '#333';

interface RegionSpan {
  name: string;
  color: string;
  startIdx: number;
  endIdx: number;
  /** Scheme position name for start, e.g. "27" in IMGT */
  startPos: string;
  /** Scheme position name for end */
  endPos: string;
  /** Source numbering scheme, e.g. "IMGT" */
  sourceScheme?: string;
  description?: string;
}

export class AnnotationTrack extends MSAHeaderTrack {
  private regions: RegionSpan[] = [];
  private posList: string[] = [];
  private _cacheVersion: number = -1;

  constructor(
    private tableCol: DG.Column<string>,
    title: string = 'Annotations',
  ) {
    super(ANNOTATION_TRACK_HEIGHT, ANNOTATION_TRACK_MIN_HEIGHT, title);
    this._rebuildRegions();
  }

  /** Returns true if there are region annotations to display. */
  hasRegions(): boolean {
    this._ensureUpToDate();
    return this.regions.length > 0;
  }

  /** Returns the region covering the given 0-based position index, or null. */
  getRegionAtPosition(posIdx: number): RegionSpan | null {
    this._ensureUpToDate();
    for (const region of this.regions) {
      if (posIdx >= region.startIdx && posIdx <= region.endIdx)
        return region;
    }
    return null;
  }

  /** Override: provide tooltip content for the annotation track. */
  public getTooltipContent(position: number, _monomer: string | null): HTMLElement | null {
    this._ensureUpToDate();
    const region = this.getRegionAtPosition(position);
    if (!region) return null;

    const scheme = region.sourceScheme ?? '';
    const posName = position < this.posList.length ? this.posList[position] : '';

    const lines: HTMLElement[] = [];
    const nameEl = document.createElement('b');
    nameEl.textContent = region.name;
    lines.push(ui.div([nameEl], {style: {fontSize: '13px', marginBottom: '4px'}}));
    if (scheme)
      lines.push(ui.divText(`Scheme: ${scheme}`, {style: {fontSize: '12px', color: '#555'}}));
    lines.push(ui.divText(`Range: ${region.startPos} \u2013 ${region.endPos}`, {style: {fontSize: '12px', color: '#555'}}));
    if (posName)
      lines.push(ui.divText(`Position: ${posName}`, {style: {fontSize: '12px', color: '#555'}}));
    if (region.description)
      lines.push(ui.divText(region.description, {style: {fontSize: '11px', color: '#888', marginTop: '2px'}}));

    return ui.divV(lines, {style: {padding: '4px'}});
  }

  draw(
    x: number, y: number, width: number, height: number,
    windowStart: number, positionWidth: number, totalPositions: number, _currentPosition: number,
  ): void {
    if (!this.ctx || !this.visible) return;
    this._ensureUpToDate();
    if (this.regions.length === 0) return;

    const g = this.ctx;
    g.save();

    // Clip to track bounds
    g.beginPath();
    g.rect(x, y, width, height);
    g.clip();

    const visibleStart = windowStart - 1; // windowStart is 1-based
    const visibleEnd = visibleStart + Math.ceil(width / positionWidth) + 2;

    for (const region of this.regions) {
      // Skip regions entirely outside the viewport
      if (region.endIdx < visibleStart || region.startIdx > visibleEnd) continue;

      const drawStart = Math.max(region.startIdx, visibleStart);
      const drawEnd = Math.min(region.endIdx, visibleEnd);

      const rx = x + (drawStart - visibleStart) * positionWidth;
      const rw = (drawEnd - drawStart + 1) * positionWidth;

      // Draw background
      g.fillStyle = region.color;
      g.globalAlpha = 0.35;
      g.fillRect(rx, y, rw, height);
      g.globalAlpha = 1.0;

      // Draw border
      g.strokeStyle = region.color;
      g.globalAlpha = 0.6;
      g.lineWidth = 1;
      g.strokeRect(rx + 0.5, y + 0.5, rw - 1, height - 1);
      g.globalAlpha = 1.0;

      // Draw label if there's enough room
      g.font = ANNOTATION_FONT;
      const labelWidth = g.measureText(region.name).width;
      if (rw > labelWidth + 4) {
        g.fillStyle = ANNOTATION_LABEL_COLOR;
        g.textAlign = 'center';
        g.textBaseline = 'middle';
        g.fillText(region.name, rx + rw / 2, y + height / 2);
      }
    }

    g.restore();
  }

  private _ensureUpToDate(): void {
    if (this._cacheVersion !== this.tableCol.version) {
      this._rebuildRegions();
      this._cacheVersion = this.tableCol.version;
    }
  }

  private _rebuildRegions(): void {
    this.regions = [];
    const annotTag = this.tableCol.getTag(bioTAGS.annotations);
    if (!annotTag) return;

    let annotations: SeqAnnotation[];
    try {
      annotations = JSON.parse(annotTag);
    } catch { return; }

    this.posList = this._getPosList();
    if (this.posList.length === 0) return;

    for (const annot of annotations) {
      if (annot.category !== AnnotationCategory.Structure) continue;
      if (annot.start == null || annot.end == null) continue;

      const startIdx = this.posList.indexOf(annot.start);
      const endIdx = this.posList.indexOf(annot.end);
      if (startIdx < 0 || endIdx < 0) continue;

      this.regions.push({
        name: annot.name,
        color: annot.color ?? this._defaultColor(annot),
        startIdx,
        endIdx,
        startPos: annot.start,
        endPos: annot.end,
        sourceScheme: annot.sourceScheme,
        description: annot.description,
      });
    }
  }

  private _getPosList(): string[] {
    const tag = this.tableCol.getTag(bioTAGS.positionNames);
    if (!tag) return [];
    return tag.split(', ');
  }

  private _defaultColor(annot: SeqAnnotation): string {
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
