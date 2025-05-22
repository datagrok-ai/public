/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
/**
 * MSAHeader.ts - Interactive MSA Header Component
 * This module provides functionality for an interactive MSA (Multiple Sequence Alignment) header
 * with position markers, slider navigation, and position selection.
 */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

// Layout Constants
const LAYOUT_CONSTANTS = {
  TITLE_HEIGHT: 12,
  TRACK_GAP: 4,
  DOTTED_CELL_HEIGHT: 30,
  SLIDER_HEIGHT: 8,
  TOP_PADDING: 5,
  DEFAULT_TRACK_HEIGHT: 45,
  MIN_TRACK_HEIGHT: 35
} as const;

// WebLogo Constants
const WEBLOGO_CONSTANTS = {
  PADDING: 2,
  MIN_LETTER_HEIGHT: 4,
  SEPARATOR_WIDTH: 1,
  MIN_GRID_COLUMNS: 45
} as const;

// Typography Constants
const FONTS = {
  TITLE: 'bold 10px sans-serif',
  POSITION_LABELS: '12px monospace',
  CONSERVATION_TEXT: '9px monospace',
  TOOLTIP_MAIN: '13px',
  TOOLTIP_SECONDARY: '12px',
  TOOLTIP_SMALL: '11px'
} as const;

// Color Constants
const COLORS = {
  TITLE_TEXT: '#333333',
  SELECTION_HIGHLIGHT: 'rgba(60, 177, 115, 0.1)',
  SELECTION_STRONG: 'rgba(60, 177, 115, 0.2)',
  SELECTION_CONNECTION: 'rgba(60, 177, 115, 0.4)',
  BACKGROUND_LIGHT: 'rgba(240, 240, 240, 0.5)',
  BORDER_LIGHT: 'rgba(100, 100, 100, 0.3)',
  SEPARATOR_LIGHT: 'rgba(255, 255, 255, 0.4)',
  HOVER_SHADOW: 'rgba(255, 255, 255, 0.8)',
  SLIDER_DEFAULT: 'rgba(220, 220, 220, 0.4)',
  SLIDER_WINDOW: 'rgba(150, 150, 150, 0.5)',
  POSITION_DOT: '#999999',
  SLIDER_MARKER: '#3CB173'
} as const;

// Tooltip Style Template
const TOOLTIP_STYLE = {
  color: '#4A4A49',
  background: '#FFFBCC',
  borderRadius: '3px',
  padding: '8px',
  fontSize: '12px',
  fontFamily: 'sans-serif',
  minWidth: '120px'
} as const;

// Position range interface
interface WindowRange {
  start: number;
  end: number;
}

// Configuration options interface
interface MSAHeaderOptions {
  canvas: HTMLCanvasElement;
  totalPositions?: number;
  positionWidth?: number;
  headerHeight?: number;
  sliderHeight?: number;
  currentPosition?: number;
  windowStartPosition?: number;
  cellBackground?: boolean;
  sliderColor?: string;
  x?: number;
  y?: number;
  width?: number;
  height?: number;
  onPositionChange?: (position: number, range: WindowRange) => void;
}

interface Preventable {
  preventDefault: () => void;
}

// Internal state interface
interface MSAHeaderState {
  isDragging: boolean;
  dragStartX: number;
}

/**
 * Base class for all MSA header tracks
 */
export abstract class MSAHeaderTrack {
  protected ctx: CanvasRenderingContext2D | null = null;
  protected visible: boolean = true;
  protected height: number;
  protected minHeight: number;
  protected defaultHeight: number;
  protected title: string = '';
  protected tooltipEnabled: boolean = false;
  protected tooltipContent: ((position: number, data: any) => HTMLElement) | null = null;

  constructor(height: number = LAYOUT_CONSTANTS.DEFAULT_TRACK_HEIGHT,
    minHeight: number = LAYOUT_CONSTANTS.MIN_TRACK_HEIGHT,
    title: string = '') {
    this.height = height;
    this.defaultHeight = height;
    this.minHeight = minHeight + (title ? LAYOUT_CONSTANTS.TITLE_HEIGHT : 0);
    this.title = title;
  }

  /**
   * Initialize the track with a canvas context
   */
  init(ctx: CanvasRenderingContext2D): void {
    this.ctx = ctx;
  }

  public getMonomerAt(x: number, y: number, position: number): string | null {
    return null;
  }

  public enableTooltip(enabled: boolean): void {
    this.tooltipEnabled = enabled;
  }

  public setTooltipContentGenerator(contentGenerator: (position: number, data: any) => HTMLElement): void {
    this.tooltipContent = contentGenerator;
  }

  public getTooltipContent(position: number): HTMLElement | null {
    if (!this.tooltipEnabled || !this.tooltipContent) return null;
    const positionData = this.getPositionData(position);
    return this.tooltipContent(position, positionData);
  }

  /**
   * Get data for a specific position (to be implemented by subclasses)
   */
  protected getPositionData(position: number): any {
    return null;
  }

  setVisible(visible: boolean): void {
    this.visible = visible;
  }

  getHeight(): number {
    return this.visible ? this.height : 0;
  }

  getDefaultHeight(): number {
    return this.defaultHeight;
  }

  getMinHeight(): number {
    return this.minHeight;
  }

  setHeight(height: number): void {
    this.height = Math.max(this.minHeight, height);
  }

  resetHeight(): void {
    this.height = this.defaultHeight;
  }

  isVisible(): boolean {
    return this.visible;
  }

  setTitle(title: string): void {
    this.title = title;
  }

  getTitle(): string {
    return this.title;
  }

  /**
   * Draw track title
   */
  protected drawTitle(x: number, y: number, width: number): void {
    if (!this.ctx || !this.title) return;

    this.ctx.fillStyle = COLORS.TITLE_TEXT;
    this.ctx.font = FONTS.TITLE;
    this.ctx.textAlign = 'left';
    this.ctx.textBaseline = 'top';
    this.ctx.fillText(this.title, x + 5, y + 2);
  }

  /**
   * Draw the track at the specified position
   */
  abstract draw(
    x: number, y: number, width: number, height: number,
    windowStart: number, positionWidth: number, totalPositions: number, currentPosition: number
  ): void;
}

export class WebLogoTrack extends MSAHeaderTrack {
  private data: Map<number, Map<string, number>> = new Map();
  private monomerLib: any = null;
  private biotype: string = 'PEPTIDE';
  private hoveredPosition: number = -1;
  private hoveredMonomer: string | null = null;

  // Fallback color scheme if monomerLib doesn't work
  private readonly aminoAcidColors: Record<string, string> = {
    'A': '#80a0f0', 'R': '#f01505', 'N': '#00ff00', 'D': '#c048c0', 'C': '#f08080',
    'Q': '#00ff00', 'E': '#c048c0', 'G': '#f09048', 'H': '#15a4a4', 'I': '#80a0f0',
    'L': '#80a0f0', 'K': '#f01505', 'M': '#80a0f0', 'F': '#80a0f0', 'P': '#ffff00',
    'S': '#00ff00', 'T': '#00ff00', 'W': '#80a0f0', 'Y': '#15a4a4', 'V': '#80a0f0'
  };

  constructor(data: Map<number, Map<string, number>> = new Map(),
    height: number = LAYOUT_CONSTANTS.DEFAULT_TRACK_HEIGHT,
    _colorScheme: string = '',
    title: string = 'WebLogo') {
    super(height, LAYOUT_CONSTANTS.DEFAULT_TRACK_HEIGHT, title);
    this.data = data;
    this.visible = data.size > 0;
  }

  public setHovered(position: number, monomer: string | null): void {
    this.hoveredPosition = position;
    this.hoveredMonomer = monomer;
  }

  protected getPositionData(position: number): Map<string, number> | null {
    return this.data.get(position) || null;
  }

  public setupDefaultTooltip(): void {
    this.enableTooltip(true);
    this.setTooltipContentGenerator((position: number, data: Map<string, number>) => {
      return this.createTooltipContent(position, data);
    });
  }

  private createTooltipContent(position: number, data: Map<string, number>): HTMLElement {
    const container = ui.div([], {style: {...TOOLTIP_STYLE}});

    // Position header
    const positionDiv = ui.div([ui.divText(`Position: ${position + 1}`)], {
      style: {fontWeight: 'bold', marginBottom: '6px', fontSize: FONTS.TOOLTIP_MAIN}
    });
    container.appendChild(positionDiv);

    if (data && data.size > 0) {
      const gridContainer = this.createFrequencyGrid(data);
      container.appendChild(gridContainer);
    } else {
      container.appendChild(ui.divText('No data available', {
        style: {fontStyle: 'italic', color: '#666'}
      }));
    }

    return container;
  }

  private createFrequencyGrid(data: Map<string, number>): HTMLElement {
    const sortedResidues = Array.from(data.entries()).sort((a, b) => b[1] - a[1]);

    const gridContainer = ui.div([], {
      style: {
        display: 'grid',
        gridTemplateColumns: `repeat(auto-fit, minmax(${WEBLOGO_CONSTANTS.MIN_GRID_COLUMNS}px, 1fr))`,
        gap: '3px',
        marginTop: '4px'
      }
    });

    for (const [residue, freq] of sortedResidues) {
      const cellDiv = this.createFrequencyCell(residue, freq);
      gridContainer.appendChild(cellDiv);
    }

    return gridContainer;
  }

  private createFrequencyCell(residue: string, freq: number): HTMLElement {
    const cellDiv = ui.div([], {
      style: {
        backgroundColor: this.getMonomerTextColorForBackground(residue),
        color: this.getContrastColor(this.getMonomerTextColorForBackground(residue)),
        textAlign: 'center',
        padding: '4px 2px',
        borderRadius: '2px',
        fontSize: FONTS.TOOLTIP_SMALL,
        fontWeight: 'bold',
        border: '1px solid rgba(0,0,0,0.1)'
      }
    });

    const residueText = ui.div([ui.divText(residue)], {
      style: {fontSize: FONTS.TOOLTIP_SECONDARY, fontWeight: 'bold', lineHeight: '1'}
    });
    const freqText = ui.div([ui.divText(`${(freq * 100).toFixed(0)}%`)], {
      style: {fontSize: '10px', lineHeight: '1', marginTop: '1px'}
    });

    cellDiv.appendChild(residueText);
    cellDiv.appendChild(freqText);
    return cellDiv;
  }

  setMonomerLib(monomerLib: any): void {
    this.monomerLib = monomerLib;
  }

  setBiotype(biotype: string): void {
    this.biotype = biotype;
  }

  /**
   * Calculate which monomer is at the specified coordinates within a WebLogo column
   */
  public getMonomerAt(x: number, y: number, position: number): string | null {
    if (!this.ctx || !this.visible || this.data.size === 0) return null;

    const residueFreqs = this.data.get(position);
    if (!residueFreqs || residueFreqs.size === 0) return null;

    const titleHeight = LAYOUT_CONSTANTS.TITLE_HEIGHT + 2;
    if (y < titleHeight) return null;

    const sortedResidues = Array.from(residueFreqs.entries()).sort((a, b) => b[1] - a[1]);
    const relativeY = y - titleHeight;
    const totalFreq = sortedResidues.reduce((sum, [_, freq]) => sum + freq, 0);
    const columnHeight = this.height - titleHeight;

    let currentY = 0;
    for (const [residue, freq] of sortedResidues) {
      const letterHeight = freq * columnHeight / totalFreq;
      if (relativeY >= currentY && relativeY < currentY + letterHeight)
        return residue;

      currentY += letterHeight;
    }

    return null;
  }

  updateData(data: Map<number, Map<string, number>>): void {
    this.data = data;
    this.visible = data.size > 0;
  }

  /**
   * Draw the WebLogo track
   */
  draw(x: number, y: number, width: number, height: number, windowStart: number,
    positionWidth: number, totalPositions: number, currentPosition: number): void {
    if (!this.ctx || !this.visible || this.data.size === 0) return;

    this.drawTitle(x, y, width);

    const visiblePositionsN = Math.floor(width / positionWidth);
    const logoHeight = height - WEBLOGO_CONSTANTS.PADDING * 2;
    const effectiveLogoHeight = logoHeight - LAYOUT_CONSTANTS.TITLE_HEIGHT - 2;
    const columnY = y + LAYOUT_CONSTANTS.TITLE_HEIGHT + WEBLOGO_CONSTANTS.PADDING;

    for (let i = 0; i < visiblePositionsN; i++) {
      this.drawWebLogoColumn(i, x, y, width, height, windowStart, positionWidth,
        totalPositions, currentPosition, columnY, effectiveLogoHeight);
    }
  }

  private drawWebLogoColumn(i: number, x: number, y: number, width: number, height: number,
    windowStart: number, positionWidth: number, totalPositions: number,
    currentPosition: number, columnY: number, effectiveLogoHeight: number): void {
    const position = windowStart + i - 1;
    if (position < 0 || position >= totalPositions) return;

    const residueFreqs = this.data.get(position);
    if (!residueFreqs || residueFreqs.size === 0) return;

    const posX = x + (i * positionWidth);
    const cellWidth = positionWidth - WEBLOGO_CONSTANTS.PADDING;

    // Highlight selected position
    if (windowStart + i === currentPosition) {
      this.ctx!.fillStyle = COLORS.SELECTION_HIGHLIGHT;
      this.ctx!.fillRect(posX, y, positionWidth, height);
    }

    // Draw column background
    this.ctx!.fillStyle = COLORS.BACKGROUND_LIGHT;
    this.ctx!.fillRect(posX + WEBLOGO_CONSTANTS.PADDING / 2, columnY, cellWidth, effectiveLogoHeight);

    this.drawLettersInColumn(position, posX, cellWidth, columnY, effectiveLogoHeight, residueFreqs);
    this.drawColumnBorder(posX, columnY, cellWidth, effectiveLogoHeight);
  }

  private drawLettersInColumn(position: number, posX: number, cellWidth: number,
    columnY: number, effectiveLogoHeight: number,
    residueFreqs: Map<string, number>): void {
    const sortedResidues = Array.from(residueFreqs.entries()).sort((a, b) => b[1] - a[1]);
    const totalFreq = sortedResidues.reduce((sum, [_, freq]) => sum + freq, 0);
    const scaleFactor = Math.min(1, effectiveLogoHeight / (totalFreq * effectiveLogoHeight));

    let currentY = columnY;
    const columnBottom = columnY + effectiveLogoHeight;

    for (const [residue, freq] of sortedResidues) {
      const rawLetterHeight = freq * effectiveLogoHeight * scaleFactor;
      const letterHeight = Math.max(WEBLOGO_CONSTANTS.MIN_LETTER_HEIGHT, Math.floor(rawLetterHeight));

      if (letterHeight < WEBLOGO_CONSTANTS.MIN_LETTER_HEIGHT) continue;

      const isHovered = position === this.hoveredPosition && residue === this.hoveredMonomer;
      const finalHeight = Math.min(letterHeight, columnBottom - currentY);

      if (finalHeight < WEBLOGO_CONSTANTS.MIN_LETTER_HEIGHT) break;

      this.drawLetter(residue, posX + WEBLOGO_CONSTANTS.PADDING / 2, currentY, cellWidth, finalHeight, isHovered);
      currentY += finalHeight;
    }
  }

  private drawColumnBorder(posX: number, columnY: number, cellWidth: number, effectiveLogoHeight: number): void {
    this.ctx!.strokeStyle = COLORS.BORDER_LIGHT;
    this.ctx!.lineWidth = 1;
    this.ctx!.strokeRect(posX + WEBLOGO_CONSTANTS.PADDING / 2, columnY, cellWidth, effectiveLogoHeight);
  }

  private drawLetter(letter: string, x: number, y: number, width: number, height: number, isHovered: boolean = false): void {
    if (!this.ctx) return;

    const backgroundColor = this.getMonomerTextColorForBackground(letter);
    const textColor = this.getContrastColor(backgroundColor);

    // Fill background
    this.ctx.fillStyle = backgroundColor;
    this.ctx.fillRect(x, y, width, height);

    // Apply hover effect
    if (isHovered) {
      this.ctx.shadowColor = COLORS.HOVER_SHADOW;
      this.ctx.shadowBlur = 8;
      this.ctx.strokeStyle = 'white';
      this.ctx.lineWidth = 2;
      this.ctx.strokeRect(x, y, width, height);
      this.ctx.shadowBlur = 0;
    }

    this.drawLetterSeparators(x, y, width, height);
    this.drawLetterText(letter, x, y, width, height, textColor);
  }

  private drawLetterSeparators(x: number, y: number, width: number, height: number): void {
    this.ctx!.strokeStyle = COLORS.SEPARATOR_LIGHT;
    this.ctx!.lineWidth = WEBLOGO_CONSTANTS.SEPARATOR_WIDTH;

    if (y + height < y + this.ctx!.canvas.height) {
      this.ctx!.beginPath();
      this.ctx!.moveTo(x, y + height);
      this.ctx!.lineTo(x + width, y + height);
      this.ctx!.stroke();
    }
  }

  private drawLetterText(letter: string, x: number, y: number, width: number, height: number, textColor: string): void {
    const fontSize = Math.min(height * 0.8, width * 0.8);
    if (fontSize >= 7) {
      this.ctx!.fillStyle = textColor;
      this.ctx!.font = `bold ${fontSize}px sans-serif`;
      this.ctx!.textAlign = 'center';
      this.ctx!.textBaseline = 'middle';
      this.ctx!.fillText(letter, x + width / 2, y + height / 2);
    }
  }

  private getMonomerTextColorForBackground(letter: string): string {
    if (this.monomerLib && typeof this.monomerLib.getMonomerTextColor === 'function') {
      try {
        const textColor = this.monomerLib.getMonomerTextColor(this.biotype, letter);
        if (textColor) return textColor;
      } catch (e) {
        console.warn('Error getting text color from monomerLib:', e);
      }
    }

    if (this.monomerLib && typeof this.monomerLib.getMonomerColor === 'function') {
      try {
        const color = this.monomerLib.getMonomerColor(this.biotype, letter);
        if (color) return color;
      } catch (e) {
        console.warn('Error getting color from monomerLib:', e);
      }
    }

    return this.aminoAcidColors[letter] || '#CCCCCC';
  }

  private getContrastColor(color: string): string {
    try {
      let r; let g; let b;
      if (color.startsWith('rgb')) {
        const rgbMatch = color.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*[\d.]+)?\)/);
        if (rgbMatch) {
          r = parseInt(rgbMatch[1], 10);
          g = parseInt(rgbMatch[2], 10);
          b = parseInt(rgbMatch[3], 10);
        } else return '#000000';
      } else if (color.startsWith('#')) {
        const hex = color.slice(1);
        r = parseInt(hex.slice(0, 2), 16);
        g = parseInt(hex.slice(2, 4), 16);
        b = parseInt(hex.slice(4, 6), 16);
      } else return '#000000';

      const yiq = ((r * 299) + (g * 587) + (b * 114)) / 1000;
      return (yiq >= 128) ? '#000000' : '#FFFFFF';
    } catch (e) {
      return '#000000';
    }
  }
}

/**
 * Track for displaying conservation bars
 */
export class ConservationTrack extends MSAHeaderTrack {
  private data: number[];
  private colorScheme: 'default' | 'rainbow' | 'heatmap';

  constructor(data: number[],
    height: number = LAYOUT_CONSTANTS.DEFAULT_TRACK_HEIGHT,
    colorScheme: 'default' | 'rainbow' | 'heatmap' = 'default',
    title: string = 'Conservation') {
    super(height, LAYOUT_CONSTANTS.MIN_TRACK_HEIGHT, title);
    this.data = data;
    this.colorScheme = colorScheme;
    this.visible = data.length > 0;
  }

  updateData(data: number[]): void {
    this.data = data;
    this.visible = data.length > 0;
  }

  /**
   * Draw the conservation track
   */
  draw(x: number, y: number, width: number, height: number, windowStart: number,
    positionWidth: number, totalPositions: number, currentPosition: number): void {
    if (!this.ctx || !this.visible || this.data.length === 0) return;

    this.drawTitle(x, y, width);

    const barplotTop = y + LAYOUT_CONSTANTS.TITLE_HEIGHT;
    const barplotHeight = height - LAYOUT_CONSTANTS.TITLE_HEIGHT;
    const visiblePositionsN = Math.floor(width / positionWidth);

    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > totalPositions) break;

      const posX = x + (i * positionWidth);
      const cellWidth = positionWidth;
      const cellCenterX = posX + cellWidth / 2;

      if (position - 1 < this.data.length) {
        this.drawConservationBar(position - 1, posX, cellWidth, cellCenterX, barplotTop, barplotHeight);

        if (position === currentPosition) {
          this.ctx.fillStyle = COLORS.SELECTION_HIGHLIGHT;
          this.ctx.fillRect(posX, barplotTop, cellWidth, barplotHeight);
        }
      }
    }
  }

  private drawConservationBar(posIndex: number, x: number, cellWidth: number, cellCenterX: number,
    barplotTop: number, barplotHeight: number): void {
    if (!this.ctx) return;

    const conservation = this.data[posIndex];
    const barplotPadding = WEBLOGO_CONSTANTS.PADDING;

    // Draw bar background
    this.ctx.fillStyle = COLORS.BACKGROUND_LIGHT;
    this.ctx.fillRect(x + barplotPadding, barplotTop, cellWidth - barplotPadding * 2, barplotHeight);

    // Determine bar color based on color scheme
    let barColor = '#3CB173';
    if (this.colorScheme === 'default') {
      if (conservation < 0.5) barColor = '#E74C3C';
      else if (conservation < 0.75) barColor = '#F39C12';
    } else if (this.colorScheme === 'rainbow') {
      if (conservation < 0.2) barColor = '#E74C3C';
      else if (conservation < 0.4) barColor = '#FF7F00';
      else if (conservation < 0.6) barColor = '#FFFF00';
      else if (conservation < 0.8) barColor = '#00FF00';
      else barColor = '#0000FF';
    } else if (this.colorScheme === 'heatmap') {
      const intensity = Math.round(conservation * 255);
      barColor = `rgb(255, ${intensity}, ${intensity})`;
    }

    const barHeight = conservation * barplotHeight;
    this.ctx.fillStyle = barColor;
    this.ctx.fillRect(x + barplotPadding, barplotTop + barplotHeight - barHeight,
      cellWidth - barplotPadding * 2, barHeight);

    // Add outline
    this.ctx.strokeStyle = COLORS.BORDER_LIGHT;
    this.ctx.lineWidth = 1;
    this.ctx.strokeRect(x + barplotPadding, barplotTop, cellWidth - barplotPadding * 2, barplotHeight);

    // Add conservation percentage text if cell is wide enough
    if (cellWidth > 20) {
      this.ctx.fillStyle = COLORS.TITLE_TEXT;
      this.ctx.font = FONTS.CONSERVATION_TEXT;
      this.ctx.textAlign = 'center';
      this.ctx.textBaseline = 'middle';
      this.ctx.fillText(`${Math.round(conservation * 100)}%`, cellCenterX, barplotTop + barplotHeight / 2);
    }
  }
}

/**
 * Main MSA header class that manages dotted cells and tracks
 */
export class MSAScrollingHeader {
  private config: Required<MSAHeaderOptions>;
  private state: MSAHeaderState;
  private canvas: HTMLCanvasElement | null = null;
  private ctx: CanvasRenderingContext2D | null = null;
  private eventElement: HTMLDivElement;

  private tracks: Map<string, MSAHeaderTrack> = new Map();

  // Tooltip and hover state
  private currentHoverPosition: number = -1;
  private currentHoverTrack: string | null = null;
  private previousHoverPosition: number = -1;
  private previousHoverTrack: string | null = null;
  private previousHoverMonomer: string | null = null;

  // Selection state
  private dataFrame: DG.DataFrame | null = null;
  private seqHandler: any = null;
  private seqColumn: DG.Column<string> | null = null;
  private onSelectionCallback: ((position: number, monomer: string) => void) | null = null;

  constructor(options: MSAHeaderOptions) {
    this.config = {
      x: options.x || 0,
      y: options.y || 0,
      width: options.width || 0,
      height: options.height || 0,
      windowStartPosition: options.windowStartPosition || 1,
      positionWidth: options.positionWidth || 15,
      totalPositions: options.totalPositions || 5000,
      headerHeight: options.headerHeight || 50,
      sliderHeight: options.sliderHeight || LAYOUT_CONSTANTS.SLIDER_HEIGHT,
      currentPosition: options.currentPosition || 1,
      cellBackground: options.cellBackground !== undefined ? options.cellBackground : true,
      sliderColor: options.sliderColor || COLORS.SLIDER_DEFAULT,
      onPositionChange: options.onPositionChange || ((_, __) => { }),
      ...options
    };

    this.eventElement = ui.div();
    this.eventElement.style.position = 'absolute';
    this.config.canvas.parentElement?.appendChild(this.eventElement);

    this.state = {isDragging: false, dragStartX: 0};
    this.setupEventListeners();
    this.init();
  }

  /**
   * Clear hover states for all tracks. Only triggers redraw if there was an active hover state.
   */
  private clearHoverStates(): void {
    const hasHoverState = this.previousHoverPosition !== -1 ||
                         this.previousHoverTrack !== null ||
                         this.previousHoverMonomer !== null;

    if (hasHoverState) {
      this.previousHoverPosition = -1;
      this.previousHoverTrack = null;
      this.previousHoverMonomer = null;

      this.tracks.forEach((track) => {
        if (track instanceof WebLogoTrack)
          track.setHovered(-1, null);
      });

      window.requestAnimationFrame(() => this.draw(
        this.config.x, this.config.y, this.config.width, this.config.height,
        this.config.currentPosition, this.config.windowStartPosition, {preventDefault: () => {}}
      ));
    }
  }

  public setSelectionData(dataFrame: DG.DataFrame, seqColumn: DG.Column<string>, seqHandler: any,
    callback?: (position: number, monomer: string) => void): void {
    this.dataFrame = dataFrame;
    this.seqColumn = seqColumn;
    this.seqHandler = seqHandler;
    this.onSelectionCallback = callback || null;
  }

  public setupTooltipHandling(): void {
    this.eventElement.addEventListener('mousemove', this.handleTooltipMouseMove.bind(this));
    this.eventElement.addEventListener('mouseleave', this.handleTooltipMouseLeave.bind(this));
  }

  /**
   * Handle mouse movement for tooltips and hover effects.
   * Only triggers redraws when hover state actually changes to prevent flickering.
   */
  private handleTooltipMouseMove(e: MouseEvent): void {
    if (!this.isValid) return;

    const {x, y} = this.getCoords(e);
    const cellWidth = this.config.positionWidth;
    const hoveredCellIndex = Math.floor(x / cellWidth);
    const windowStart = this.config.windowStartPosition;
    const hoveredPosition = windowStart + hoveredCellIndex - 1;

    if (hoveredPosition < 0 || hoveredPosition >= this.config.totalPositions) {
      this.hideTooltip();
      this.clearHoverStates();
      return;
    }

    const sliderTop = this.config.headerHeight - LAYOUT_CONSTANTS.SLIDER_HEIGHT;
    const dottedCellsTop = sliderTop - LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT;

    if (y >= dottedCellsTop) {
      this.hideTooltip();
      this.clearHoverStates();
      return;
    }

    // Find which track is being hovered
    let trackY = 0;
    let hoveredTrackId: string | null = null;
    let trackRelativeY = 0;

    for (const [id, track] of this.tracks.entries()) {
      if (!track.isVisible()) continue;

      const trackHeight = track.getHeight();
      if (y >= trackY && y < trackY + trackHeight) {
        hoveredTrackId = id;
        trackRelativeY = y - trackY;
        break;
      }

      trackY += trackHeight + LAYOUT_CONSTANTS.TRACK_GAP;
    }

    // Check for hovered monomer
    let currentHoverMonomer: string | null = null;
    if (hoveredTrackId) {
      const track = this.tracks.get(hoveredTrackId);
      if (track)
        currentHoverMonomer = track.getMonomerAt(x, trackRelativeY, hoveredPosition);
    }

    // Only update and redraw if the hover state has actually changed
    const hoverStateChanged = (
      this.previousHoverPosition !== hoveredPosition ||
      this.previousHoverTrack !== hoveredTrackId ||
      this.previousHoverMonomer !== currentHoverMonomer
    );

    if (hoverStateChanged) {
      this.previousHoverPosition = hoveredPosition;
      this.previousHoverTrack = hoveredTrackId;
      this.previousHoverMonomer = currentHoverMonomer;

      if (hoveredTrackId) {
        const track = this.tracks.get(hoveredTrackId);
        if (track instanceof WebLogoTrack) {
          track.setHovered(hoveredPosition, currentHoverMonomer);

          window.requestAnimationFrame(() => this.draw(
            this.config.x, this.config.y, this.config.width, this.config.height,
            this.config.currentPosition, this.config.windowStartPosition, {preventDefault: () => {}}
          ));
        }
      }
    }

    // Handle tooltip updates
    if (hoveredPosition !== this.currentHoverPosition || hoveredTrackId !== this.currentHoverTrack) {
      this.currentHoverPosition = hoveredPosition;
      this.currentHoverTrack = hoveredTrackId;

      if (hoveredTrackId) {
        const track = this.tracks.get(hoveredTrackId);
        if (track) {
          const tooltipContent = track.getTooltipContent(hoveredPosition);
          if (tooltipContent) {
            ui.tooltip.show(tooltipContent, e.clientX + 16, e.clientY + 16);
            return;
          }
        }
      }
    }

    if (!hoveredTrackId || !currentHoverMonomer) {
      this.hideTooltip();
      this.clearHoverStates();
    }
  }

  private handleTooltipMouseLeave(): void {
    this.hideTooltip();
    this.clearHoverStates();
  }

  private hideTooltip(): void {
    this.currentHoverPosition = -1;
    this.currentHoverTrack = null;
    ui.tooltip.hide();
  }

  /**
   * Handle clicks for monomer selection
   */
  private handleSelectionClick(e: MouseEvent): void {
    if (!this.isValid || !this.dataFrame || !this.seqColumn || !this.seqHandler) return;

    const {x, y} = this.getCoords(e);
    const cellWidth = this.config.positionWidth;
    const clickedCellIndex = Math.floor(x / cellWidth);
    const windowStart = this.config.windowStartPosition;
    const clickedPosition = windowStart + clickedCellIndex - 1;

    if (clickedPosition < 0 || clickedPosition >= this.config.totalPositions) return;

    const sliderTop = this.config.headerHeight - LAYOUT_CONSTANTS.SLIDER_HEIGHT;
    const dottedCellsTop = sliderTop - LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT;

    if (y >= dottedCellsTop) return;

    // Find which track was clicked
    let trackY = 0;
    let clickedTrackId: string | null = null;
    let trackRelativeY = 0;

    for (const [id, track] of this.tracks.entries()) {
      if (!track.isVisible()) continue;

      const trackHeight = track.getHeight();
      if (y >= trackY && y < trackY + trackHeight) {
        clickedTrackId = id;
        trackRelativeY = y - trackY;
        break;
      }

      trackY += trackHeight + LAYOUT_CONSTANTS.TRACK_GAP;
    }

    if (clickedTrackId) {
      const track = this.tracks.get(clickedTrackId);
      if (track) {
        const monomer = track.getMonomerAt(x, trackRelativeY, clickedPosition);

        if (monomer) {
          if (this.onSelectionCallback) {
            this.onSelectionCallback(clickedPosition, monomer);
            return;
          }

          this.selectRowsWithMonomerAtPosition(clickedPosition, monomer);
        }
      }
    }
  }

  private selectRowsWithMonomerAtPosition(position: number, monomer: string): void {
    if (!this.dataFrame || !this.seqHandler) return;

    try {
      const selBS = DG.BitSet.create(this.dataFrame.rowCount, (rowI: number) => {
        const seqSplitted = this.seqHandler.getSplitted(rowI);
        if (!seqSplitted || position >= seqSplitted.length) return false;

        const residue = seqSplitted.getCanonical(position);
        return residue === monomer;
      });

      this.dataFrame.selection.init((i) => selBS.get(i));
    } catch (error) {
      console.error('Error selecting rows:', error);
    }
  }

  private setupEventListeners(): void {
    this.eventElement.addEventListener('mousemove', (e) => {
      if (!this.isValid) return;
      if (this.isInSliderDraggableArea(e)) this.eventElement.style.cursor = 'grab';
      else if (this.isInSliderArea(e)) this.eventElement.style.cursor = 'pointer';
      else if (this.isInHeaderArea(e)) this.eventElement.style.cursor = 'pointer';
      else this.eventElement.style.cursor = 'default';
    });

    this.eventElement.addEventListener('mousedown', this.handleMouseDown.bind(this));
    this.eventElement.addEventListener('mousemove', this.handleMouseMove.bind(this));
    this.eventElement.addEventListener('mouseup', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('mouseleave', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('click', this.handleClick.bind(this));
    this.eventElement.addEventListener('wheel', this.handleMouseWheel.bind(this));
    this.eventElement.addEventListener('click', this.handleSelectionClick.bind(this));

    window.addEventListener('keydown', this.handleKeyDown.bind(this));
  }

  private init(): void {
    this.canvas = this.config.canvas;
    if (!this.canvas) {
      console.error('Canvas not found');
      return;
    }

    const context = this.canvas.getContext('2d');
    if (!context) {
      console.error('Failed to get 2D context from canvas');
      return;
    }
    this.ctx = context;

    this.tracks.forEach((track) => track.init(context));
  }

  public addTrack(id: string, track: MSAHeaderTrack): void {
    if (this.ctx) track.init(this.ctx);
    this.tracks.set(id, track);
  }

  public removeTrack(id: string): void {
    this.tracks.delete(id);
  }

  public getTrack<T extends MSAHeaderTrack>(id: string): T | undefined {
    return this.tracks.get(id) as T | undefined;
  }

  public updateTrack<T extends MSAHeaderTrack>(id: string, updater: (track: T) => void): void {
    const track = this.getTrack<T>(id);
    if (track) updater(track);
  }

  public get isValid() {
    return !!this.canvas && !!this.ctx &&
           this.config.height >= LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT + LAYOUT_CONSTANTS.SLIDER_HEIGHT;
  }

  private calculateTotalTracksHeight(): number {
    let totalHeight = 0;
    let firstTrack = true;

    this.tracks.forEach((track) => {
      if (track.isVisible()) {
        if (!firstTrack) totalHeight += LAYOUT_CONSTANTS.TRACK_GAP;
        else firstTrack = false;
        totalHeight += track.getHeight();
      }
    });

    return totalHeight;
  }

  /**
   * Determine visible tracks based on available height and distribute extra space
   */
  private determineVisibleTracks(): void {
    const availableHeight = this.config.headerHeight -
                           (LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT + LAYOUT_CONSTANTS.SLIDER_HEIGHT + LAYOUT_CONSTANTS.TRACK_GAP);

    if (availableHeight <= 0) {
      this.tracks.forEach((track) => track.setVisible(false));
      return;
    }

    const DEFAULT_HEIGHTS = {'weblogo': 45, 'conservation': 45};

    this.tracks.forEach((track, id) => {
      const defaultHeight = DEFAULT_HEIGHTS[id as keyof typeof DEFAULT_HEIGHTS] || track.getDefaultHeight();
      track.setHeight(defaultHeight);
    });

    const prioritizedTracks = Array.from(this.tracks.entries())
      .map(([id, track]) => ({id, track}))
      .sort((a, b) => {
        const priorityA = a.id === 'weblogo' ? 2 : (a.id === 'conservation' ? 1 : 0);
        const priorityB = b.id === 'weblogo' ? 2 : (b.id === 'conservation' ? 1 : 0);
        return priorityB - priorityA;
      });

    let minRequiredHeight = 0;
    let trackCount = 0;

    for (const {track} of prioritizedTracks) {
      if (trackCount > 0) minRequiredHeight += LAYOUT_CONSTANTS.TRACK_GAP;
      minRequiredHeight += track.getMinHeight();
      trackCount++;
    }

    if (minRequiredHeight <= availableHeight) {
      prioritizedTracks.forEach(({track}) => track.setVisible(true));

      // Distribute extra space to WebLogo track
      const extraSpace = availableHeight - minRequiredHeight;
      const webLogoTrack = this.getTrack('weblogo');

      if (webLogoTrack && extraSpace > 0) {
        const currentHeight = webLogoTrack.getHeight();
        webLogoTrack.setHeight(currentHeight + extraSpace);
      }
    } else {
      this.tracks.forEach((track) => track.setVisible(false));

      let remainingHeight = availableHeight;
      let isFirstTrack = true;

      for (const {id, track} of prioritizedTracks) {
        const minTrackHeight = track.getMinHeight();
        const spaceNeeded = isFirstTrack ? minTrackHeight : minTrackHeight + LAYOUT_CONSTANTS.TRACK_GAP;

        if (remainingHeight >= spaceNeeded) {
          track.setVisible(true);
          const defaultHeight = DEFAULT_HEIGHTS[id as keyof typeof DEFAULT_HEIGHTS] || minTrackHeight;
          track.setHeight(defaultHeight);

          remainingHeight -= spaceNeeded;
          isFirstTrack = false;
        } else
          track.setVisible(false);
      }
    }
  }

  public draw(x: number, y: number, w: number, h: number, currentPos: number,
    scrollerStart: number, preventable: Preventable): void {
    const heightChanged = this.config.height !== h;

    if (this.config.x != x || this.config.y != y || this.config.width != w || this.config.height != h ||
        this.config.currentPosition != currentPos || this.config.windowStartPosition != scrollerStart) {
      Object.assign(this.config, {
        x, y, width: w, height: h,
        currentPosition: currentPos,
        windowStartPosition: scrollerStart
      });
    }

    if (!this.isValid) {
      this.eventElement.style.display = 'none';
      return;
    }

    if (heightChanged)
      this.tracks.forEach((track) => track.resetHeight());


    this.ctx!.save();
    this.ctx!.clearRect(x, y, w, h);
    this.ctx!.translate(x, y);
    this.ctx!.rect(0, 0, w, h);
    this.ctx!.clip();

    this.determineVisibleTracks();

    const sliderTop = h - LAYOUT_CONSTANTS.SLIDER_HEIGHT;
    const dottedCellsTop = sliderTop - LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT;
    const totalTracksHeight = this.calculateTotalTracksHeight();

    let trackY = Math.max(0, dottedCellsTop - totalTracksHeight - LAYOUT_CONSTANTS.TRACK_GAP);
    const visibleTrackPositions: { y: number, height: number }[] = [];

    // Get visible tracks in correct order
    const visibleTracks: Array<{id: string, track: MSAHeaderTrack}> = [];
    const conservationTrack = this.getTrack<MSAHeaderTrack>('conservation');
    if (conservationTrack && conservationTrack.isVisible())
      visibleTracks.push({id: 'conservation', track: conservationTrack});


    const webLogoTrack = this.getTrack<MSAHeaderTrack>('weblogo');
    if (webLogoTrack && webLogoTrack.isVisible())
      visibleTracks.push({id: 'weblogo', track: webLogoTrack});


    // Draw tracks
    for (const {track} of visibleTracks) {
      if (track.isVisible()) {
        const trackHeight = track.getHeight();

        track.draw(0, trackY, w, trackHeight, this.config.windowStartPosition,
          this.config.positionWidth, this.config.totalPositions, this.config.currentPosition);

        visibleTrackPositions.push({y: trackY, height: trackHeight});
        trackY += trackHeight + LAYOUT_CONSTANTS.TRACK_GAP;
      }
    }

    this.drawDottedCells(0, dottedCellsTop, w, LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT, sliderTop);
    visibleTrackPositions.push({y: dottedCellsTop, height: LAYOUT_CONSTANTS.DOTTED_CELL_HEIGHT});

    // Draw connection lines for selected position
    if (this.config.currentPosition >= 1 && this.config.currentPosition <= this.config.totalPositions) {
      const cellWidth = this.config.positionWidth;
      const position = this.config.currentPosition;
      const windowStart = this.config.windowStartPosition;
      const visibleIndex = position - windowStart;

      if (visibleIndex >= 0 && visibleIndex < Math.floor(w / cellWidth)) {
        const posX = visibleIndex * cellWidth;
        const cellCenterX = posX + cellWidth / 2;

        for (let i = 0; i < visibleTrackPositions.length - 1; i++) {
          const upperTrack = visibleTrackPositions[i];
          const lowerTrack = visibleTrackPositions[i + 1];

          this.ctx!.strokeStyle = COLORS.SELECTION_CONNECTION;
          this.ctx!.lineWidth = 1;
          this.ctx!.beginPath();
          this.ctx!.moveTo(cellCenterX, upperTrack.y + upperTrack.height);
          this.ctx!.lineTo(cellCenterX, lowerTrack.y);
          this.ctx!.stroke();
        }
      }
    }

    this.ctx!.restore();
    preventable.preventDefault();
    this.setupEventElement();
  }

  /**
   * Draw the dotted cells area with position markers
   */
  private drawDottedCells(x: number, y: number, width: number, height: number, sliderTop: number): void {
    if (!this.ctx) return;

    const totalPositions = this.config.totalPositions;
    const positionWidth = this.config.positionWidth;
    const currentPosition = this.config.currentPosition;
    const windowStart = this.config.windowStartPosition;
    const visiblePositionsN = Math.floor(width / positionWidth);
    const posIndexTop = y + LAYOUT_CONSTANTS.TOP_PADDING;

    this.drawSlider(x, sliderTop, width);

    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > totalPositions) break;

      const posX = x + (i * positionWidth);
      const cellWidth = positionWidth;
      const cellCenterX = posX + cellWidth / 2;

      if (this.config.cellBackground) {
        this.ctx.fillStyle = i % 2 === 0 ? 'rgba(248, 248, 248, 0.3)' : 'rgba(242, 242, 242, 0.2)';
        this.ctx.fillRect(posX, y, cellWidth, height);

        this.ctx.strokeStyle = 'rgba(220, 220, 220, 0.7)';
        this.ctx.beginPath();
        this.ctx.moveTo(posX, y);
        this.ctx.lineTo(posX, sliderTop);
        this.ctx.stroke();
      }

      // Draw position dot
      this.ctx.fillStyle = COLORS.POSITION_DOT;
      this.ctx.beginPath();
      this.ctx.arc(cellCenterX, posIndexTop + 5, 1, 0, Math.PI * 2);
      this.ctx.fill();

      // Draw position number for every 10th position or current position
      if (position === currentPosition || ((position === 1 || position % 10 === 0) &&
          Math.abs(position - currentPosition) > 1)) {
        this.ctx.fillStyle = COLORS.TITLE_TEXT;
        this.ctx.font = FONTS.POSITION_LABELS;
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText(position.toString(), cellCenterX, posIndexTop + 15);
      }

      // Highlight current selected position
      if (position === currentPosition) {
        this.ctx.fillStyle = COLORS.SELECTION_STRONG;
        this.ctx.fillRect(posX, y, cellWidth, height);
      }
    }
  }

  /**
   * Draw the position slider
   */
  private drawSlider(x: number, sliderTop: number, width: number): void {
    if (!this.ctx) return;

    this.ctx.fillStyle = this.config.sliderColor;
    this.ctx.fillRect(x, sliderTop, width, LAYOUT_CONSTANTS.SLIDER_HEIGHT);

    const visiblePositionsN = Math.floor(width / this.config.positionWidth);
    const windowStart = this.config.windowStartPosition;
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderWidth = this.sliderWidth;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      (windowStart - 1) / totalSliderRange * (width - sliderWidth);

    const sliderLengthPX = totalSliderRange <= 0 ? width : sliderWidth;

    this.ctx.fillStyle = COLORS.SLIDER_WINDOW;
    this.ctx.fillRect(x + sliderStartPX, sliderTop, sliderLengthPX, LAYOUT_CONSTANTS.SLIDER_HEIGHT);

    if (this.config.currentPosition >= 1 && this.config.currentPosition <= this.config.totalPositions) {
      const currentPositionRatio = (this.config.currentPosition - 1) / (this.config.totalPositions - 1);
      const notchX = Math.round(currentPositionRatio * width);

      this.ctx.fillStyle = COLORS.SLIDER_MARKER;
      this.ctx.fillRect(x + notchX - 1, sliderTop - 2, 3, LAYOUT_CONSTANTS.SLIDER_HEIGHT + 4);
    }
  }

  private setupEventElement(): void {
    this.eventElement.style.display = 'block';
    this.eventElement.style.left = `${this.config.x}px`;
    this.eventElement.style.top = `${this.config.y}px`;
    this.eventElement.style.width = `${this.config.width}px`;
    this.eventElement.style.height = `${this.config.height}px`;
  }

  getCoords(e: MouseEvent) {
    const rect = this.canvas!.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;
    const y = e.clientY - rect.top - this.config.y;
    return {x, y};
  }

  isInHeaderArea(e: MouseEvent): boolean {
    const {x, y} = this.getCoords(e);
    return x >= 0 && x <= this.config.width && y >= 0 && y <= this.config.headerHeight;
  }

  get positionWidth(): number {
    return this.config.positionWidth;
  }

  public set positionWidth(value: number) {
    this.config.positionWidth = value;
  }

  isInSliderArea(e: MouseEvent): boolean {
    const {y} = this.getCoords(e);
    const sliderTop = this.config.headerHeight - LAYOUT_CONSTANTS.SLIDER_HEIGHT;
    return y > sliderTop && y < sliderTop + LAYOUT_CONSTANTS.SLIDER_HEIGHT;
  }

  get sliderWidth(): number {
    const pseudoPositionWidth = this.config.width / this.config.totalPositions;
    const w = pseudoPositionWidth * (this.config.width / this.config.positionWidth);
    return Math.max(w, 20);
  }

  isInSliderDraggableArea(e: MouseEvent): boolean {
    const {x, y} = this.getCoords(e);
    const sliderTop = this.config.headerHeight - LAYOUT_CONSTANTS.SLIDER_HEIGHT;
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = this.config.windowStartPosition;

    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      (windowStart - 1) / totalSliderRange * (this.config.width - this.sliderWidth);

    return y > sliderTop && y < sliderTop + LAYOUT_CONSTANTS.SLIDER_HEIGHT &&
           x >= sliderStartPX && x < sliderStartPX + this.sliderWidth;
  }

  private handleMouseDown(e: MouseEvent): void {
    if (!this.isValid) return;
    const {x} = this.getCoords(e);
    if (this.isInSliderDraggableArea(e)) {
      this.state.isDragging = true;
      this.state.dragStartX = x;
      this.handleSliderDrag(x);
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
    }
  }

  private handleMouseWheel(e: WheelEvent): void {
    if (!this.isValid) return;
    if (this.isInHeaderArea(e)) {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      const delta = e.shiftKey ? Math.sign(e.deltaY) : Math.sign(e.deltaX || e.deltaY);
      const scrollSpeed = e.shiftKey ? 3 : 1;
      const newStartPosition = this.config.windowStartPosition + (delta * scrollSpeed);

      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const maxStart = this.config.totalPositions - visiblePositions + 1;
      this.config.windowStartPosition = Math.max(1, Math.min(maxStart, newStartPosition));

      if (typeof this.config.onPositionChange === 'function')
        this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  private handleMouseMove(e: MouseEvent): void {
    if (!this.state.isDragging || !this.isValid) return;

    const rect = this.canvas!.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;

    this.handleSliderDrag(x);
    e.preventDefault();
    e.stopPropagation();
    e.stopImmediatePropagation();
  }

  private handleKeyDown(e: KeyboardEvent): void {
    if (!this.isValid || this.config.currentPosition < 1) return;
    if (!document.activeElement?.contains(this.eventElement) || this.eventElement.style.display !== 'block') return;

    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      const delta = e.key === 'ArrowLeft' ? -1 : 1;
      const newPosition = Math.min(Math.max(this.config.currentPosition + delta, 1), this.config.totalPositions);

      if (newPosition === this.config.currentPosition) return;
      this.config.currentPosition = newPosition;

      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const start = this.config.windowStartPosition;
      const end = start + visiblePositions - 1;

      if (newPosition < start || newPosition > end) {
        if (delta < 0) this.config.windowStartPosition = newPosition;
        else this.config.windowStartPosition = Math.max(1, newPosition - visiblePositions + 1);
      }
    } else if (e.key === 'Escape') {
      this.config.currentPosition = -2;
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
    } else return;

    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  private handleMouseUp(): void {
    this.state.isDragging = false;
  }

  private handleSliderDrag(x: number): void {
    if (!this.isValid) return;

    const sliderWidth = this.sliderWidth;
    const canvasWidth = this.config.width - sliderWidth;
    const normalizedX = Math.max(0, Math.min(this.config.width, x));
    const fittedPositions = Math.floor(this.config.width / this.config.positionWidth);
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;

    const sliderStartPx = Math.max(0, normalizedX - sliderWidth / 2);
    const windowStart = sliderStartPx / (canvasWidth) * (totalSliderRange);
    this.config.windowStartPosition = Math.max(1, Math.min(windowStart, this.config.totalPositions - fittedPositions + 1));

    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  get headerHeight(): number {
    return this.config.headerHeight;
  }

  set headerHeight(value: number) {
    this.config.headerHeight = value;
  }

  private handleClick(e: MouseEvent): void {
    if (!this.isValid) return;

    const {x, y} = this.getCoords(e);
    const sliderTop = this.config.headerHeight - LAYOUT_CONSTANTS.SLIDER_HEIGHT;

    if (y < sliderTop && y >= 0) {
      const cellWidth = this.config.positionWidth;
      const clickedCellIndex = Math.round(x / cellWidth - 0.5);
      const windowStart = this.config.windowStartPosition;
      const clickedPosition = windowStart + clickedCellIndex;

      if (clickedPosition >= 1 && clickedPosition <= this.config.totalPositions) {
        this.config.currentPosition = clickedPosition;

        if (typeof this.config.onPositionChange === 'function')
          this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
      }
    }
  }

  public getWindowRange(): WindowRange {
    return {
      start: this.config.windowStartPosition,
      end: Math.min(this.config.totalPositions,
        this.config.windowStartPosition + Math.floor(this.config.width / this.config.positionWidth))
    };
  }

  public updateConfig(newConfig: Partial<MSAHeaderOptions>): void {
    Object.assign(this.config, newConfig);
    this.config.currentPosition = Math.min(this.config.currentPosition, this.config.totalPositions);

    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  public getCurrentPosition(): number {
    return this.config.currentPosition;
  }

  public setCurrentPosition(position: number): void {
    this.config.currentPosition = Math.max(1, Math.min(this.config.totalPositions, position));

    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }
}
