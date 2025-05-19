/* eslint-disable max-len */
/**
 * MSAHeader.ts - Interactive MSA Header Component
 * This module provides functionality for an interactive MSA (Multiple Sequence Alignment) header
 * with position markers, slider navigation, and position selection.
 */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

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
  // conservationData?: number[];
  // conservationHeight?: number;
  // conservationColorScheme?: 'default' | 'rainbow' | 'heatmap';
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

  constructor(height: number = 40, minHeight: number = 30) {
    this.height = height;
    this.minHeight = minHeight;
  }

  /**
   * Initialize the track with a canvas context
   */
  init(ctx: CanvasRenderingContext2D): void {
    this.ctx = ctx;
  }

  /**
   * Set track visibility
   */
  setVisible(visible: boolean): void {
    this.visible = visible;
  }

  /**
   * Get track height
   */
  getHeight(): number {
    return this.visible ? this.height : 0;
  }

  /**
   * Get minimum required height for this track
   */
  getMinHeight(): number {
    return this.minHeight;
  }

  /**
   * Set track height
   */
  setHeight(height: number): void {
    this.height = height;
  }

  /**
   * Whether the track is visible
   */
  isVisible(): boolean {
    return this.visible;
  }

  /**
   * Draw the track at the specified position
   * @param x X position to start drawing
   * @param y Y position to start drawing
   * @param width Width of the drawing area
   * @param height Height of the drawing area
   * @param windowStart First position visible in the window
   * @param positionWidth Width of each position
   * @param totalPositions Total number of positions
   * @param currentPosition Currently selected position
   */
  abstract draw(
    x: number,
    y: number,
    width: number,
    height: number,
    windowStart: number,
    positionWidth: number,
    totalPositions: number,
    currentPosition: number
  ): void;
}

export class WebLogoTrack extends MSAHeaderTrack {
  private data: Map<number, Map<string, number>> = new Map();
  private colorScheme: 'hydrophobicity' | 'chemistry' | 'default' = 'default';
  private padding: number = 2;
  private minLetterHeight: number = 4; // Minimum letter height to be drawn

  // Color maps for different amino acid properties
  private readonly aminoAcidColors: Record<string, Record<string, string>> = {
    // Default color scheme based on Clustal
    'default': {
      'A': '#80a0f0', 'R': '#f01505', 'N': '#00ff00', 'D': '#c048c0', 'C': '#f08080',
      'Q': '#00ff00', 'E': '#c048c0', 'G': '#f09048', 'H': '#15a4a4', 'I': '#80a0f0',
      'L': '#80a0f0', 'K': '#f01505', 'M': '#80a0f0', 'F': '#80a0f0', 'P': '#ffff00',
      'S': '#00ff00', 'T': '#00ff00', 'W': '#80a0f0', 'Y': '#15a4a4', 'V': '#80a0f0'
    },
    // Chemistry-based color scheme
    'chemistry': {
      'A': '#ccff00', 'R': '#0000ff', 'N': '#00ff00', 'D': '#ff0000', 'C': '#ffff00',
      'Q': '#00ff00', 'E': '#ff0000', 'G': '#ccff00', 'H': '#0000ff', 'I': '#ccff00',
      'L': '#ccff00', 'K': '#0000ff', 'M': '#ccff00', 'F': '#ccff00', 'P': '#ffff00',
      'S': '#00ff00', 'T': '#00ff00', 'W': '#ccff00', 'Y': '#00ff00', 'V': '#ccff00'
    },
    // Hydrophobicity scale
    'hydrophobicity': {
      'A': '#ad0052', 'R': '#0000ff', 'N': '#0c00f3', 'D': '#0c00f3', 'C': '#c2003d',
      'Q': '#0c00f3', 'E': '#0c00f3', 'G': '#6a0095', 'H': '#1500ea', 'I': '#ff0000',
      'L': '#ea0015', 'K': '#0000ff', 'M': '#b0004f', 'F': '#cb0034', 'P': '#4600b9',
      'S': '#5e00a1', 'T': '#61009e', 'W': '#5e0093', 'Y': '#4f00b0', 'V': '#f60009'
    }
  };

  constructor(
    data: Map<number, Map<string, number>> = new Map(),
    height: number = 40, // Match conservation barplot height by default
    colorScheme: 'hydrophobicity' | 'chemistry' | 'default' = 'default'
  ) {
    super(height, 40); // Min height also set to 40px to match conservation track
    this.data = data;
    this.colorScheme = colorScheme;
    this.visible = data.size > 0;
  }

  /**
   * Update the WebLogo data
   */
  updateData(data: Map<number, Map<string, number>>): void {
    this.data = data;
    this.visible = data.size > 0;
  }

  /**
   * Set color scheme
   */
  setColorScheme(scheme: 'hydrophobicity' | 'chemistry' | 'default'): void {
    this.colorScheme = scheme;
  }

  /**
   * Draw the WebLogo track
   */
  draw(
    x: number,
    y: number,
    width: number,
    height: number,
    windowStart: number,
    positionWidth: number,
    totalPositions: number,
    currentPosition: number
  ): void {
    if (!this.ctx || !this.visible || this.data.size === 0) return;

    const visiblePositionsN = Math.floor(width / positionWidth);
    const logoHeight = height - this.padding * 2;
    const titleHeight = 12;
    const effectiveLogoHeight = logoHeight - titleHeight - 2;

    // Draw title
    this.ctx.fillStyle = '#333333';
    this.ctx.font = 'bold 10px sans-serif';
    this.ctx.textAlign = 'left';
    this.ctx.textBaseline = 'top';
    // this.ctx.fillText('WebLogo', x + 5, y + this.padding);

    // Draw WebLogo columns
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i - 1; // Subtract 1 to convert to 0-based index
      if (position < 0 || position >= totalPositions) continue;

      const residueFreqs = this.data.get(position);
      if (!residueFreqs || residueFreqs.size === 0) continue;

      const posX = x + (i * positionWidth);
      const cellWidth = positionWidth - this.padding;
      const columnY = y + titleHeight + this.padding;

      // Highlight current selected position
      if (windowStart + i === currentPosition) {
        this.ctx.fillStyle = 'rgba(60, 177, 115, 0.1)';
        this.ctx.fillRect(posX, y, positionWidth, height);
      }

      // Draw background for the column
      this.ctx.fillStyle = 'rgba(240, 240, 240, 0.5)';
      this.ctx.fillRect(
        posX + this.padding / 2,
        columnY,
        cellWidth,
        effectiveLogoHeight
      );

      // Sort residues by frequency (highest to lowest)
      const sortedResidues = Array.from(residueFreqs.entries())
        .sort((a, b) => b[1] - a[1]);

      // Draw each residue letter proportional to its frequency
      // But ensure all columns have the same total height
      let currentY = columnY;

      for (const [residue, freq] of sortedResidues) {
        // Calculate letter height proportional to frequency
        const letterHeight = Math.max(this.minLetterHeight, Math.floor(freq * effectiveLogoHeight));

        // Skip very small letters
        if (letterHeight < this.minLetterHeight) continue;

        // Draw the letter
        this.drawLetter(
          residue,
          posX + this.padding / 2,
          currentY,
          cellWidth,
          letterHeight
        );

        // Move current Y position down for next letter
        currentY += letterHeight;

        // Stop if we've reached the bottom of the column
        if (currentY >= columnY + effectiveLogoHeight) break;
      }

      // Draw column border
      this.ctx.strokeStyle = 'rgba(100, 100, 100, 0.3)';
      this.ctx.lineWidth = 1;
      this.ctx.strokeRect(
        posX + this.padding / 2,
        columnY,
        cellWidth,
        effectiveLogoHeight
      );
    }
  }

  /**
   * Draw a letter in the WebLogo
   */
  private drawLetter(
    letter: string,
    x: number,
    y: number,
    width: number,
    height: number
  ): void {
    if (!this.ctx) return;

    // Get color for this amino acid
    const colorMap = this.aminoAcidColors[this.colorScheme] || this.aminoAcidColors.default;
    const color = colorMap[letter] || '#888888'; // Default gray for unknown residues

    // Fill background with amino acid color
    this.ctx.fillStyle = color;
    this.ctx.fillRect(x, y, width, height);

    // Draw letter
    this.ctx.fillStyle = this.getContrastColor(color);

    // Adjust font size based on letter box size
    const fontSize = Math.min(height * 0.8, width * 0.8);
    if (fontSize >= 7) { // Only draw text if it will be readable
      this.ctx.font = `bold ${fontSize}px sans-serif`;
      this.ctx.textAlign = 'center';
      this.ctx.textBaseline = 'middle';
      this.ctx.fillText(letter, x + width / 2, y + height / 2);
    }
  }

  /**
   * Calculate a contrasting text color (black or white) based on background color
   */
  private getContrastColor(hexColor: string): string {
    // Convert hex to RGB
    const r = parseInt(hexColor.slice(1, 3), 16);
    const g = parseInt(hexColor.slice(3, 5), 16);
    const b = parseInt(hexColor.slice(5, 7), 16);

    // Calculate perceived brightness (YIQ formula)
    const yiq = ((r * 299) + (g * 587) + (b * 114)) / 1000;

    // Return black or white based on brightness
    return (yiq >= 128) ? 'black' : 'white';
  }
}
/**
 * Track for displaying conservation bars
 */
export class ConservationTrack extends MSAHeaderTrack {
  private data: number[];
  private colorScheme: 'default' | 'rainbow' | 'heatmap';
  private padding: number = 2;

  constructor(data: number[], height: number = 40, colorScheme: 'default' | 'rainbow' | 'heatmap' = 'default') {
    super(height, 30);
    this.data = data;
    this.colorScheme = colorScheme;
    this.visible = data.length > 0;
  }

  /**
   * Update the conservation data
   */
  updateData(data: number[]): void {
    this.data = data;
    this.visible = data.length > 0;
  }

  /**
   * Draw the conservation track
   */
  draw(
    x: number,
    y: number,
    width: number,
    height: number,
    windowStart: number,
    positionWidth: number,
    totalPositions: number,
    currentPosition: number
  ): void {
    if (!this.ctx || !this.visible || this.data.length === 0) return;

    const barplotTop = y;
    const barplotHeight = height;
    const visiblePositionsN = Math.floor(width / positionWidth);

    // Draw position conservation bars
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > totalPositions) break;

      const posX = x + (i * positionWidth);
      const cellWidth = positionWidth;
      const cellCenterX = posX + cellWidth / 2;

      if (position - 1 < this.data.length) {
        this.drawConservationBar(
          position - 1,
          posX,
          cellWidth,
          cellCenterX,
          barplotTop,
          barplotHeight,
          this.padding
        );

        // Highlight current selected position
        if (position === currentPosition) {
          this.ctx.fillStyle = 'rgba(60, 177, 115, 0.1)';
          this.ctx.fillRect(
            posX,
            barplotTop,
            cellWidth,
            barplotHeight
          );
        }
      }
    }
  }

  /**
   * Draw a conservation bar for a specific position
   */
  private drawConservationBar(
    posIndex: number,
    x: number,
    cellWidth: number,
    cellCenterX: number,
    barplotTop: number,
    barplotHeight: number,
    barplotPadding: number
  ): void {
    if (!this.ctx) return;

    const conservation = this.data[posIndex];

    // Draw bar background
    this.ctx.fillStyle = 'rgba(240, 240, 240, 0.5)';
    this.ctx.fillRect(
      x + barplotPadding,
      barplotTop,
      cellWidth - barplotPadding * 2,
      barplotHeight
    );

    // Determine bar color based on color scheme
    let barColor = '#3CB173'; // Default green for high conservation

    if (this.colorScheme === 'default') {
      // Default scheme: green (high), yellow (medium), red (low)
      if (conservation < 0.5) {
        barColor = '#E74C3C'; // Red for low conservation (<50%)
      } else if (conservation < 0.75) {
        barColor = '#F39C12'; // Yellow for medium conservation (50-75%)
      }
    } else if (this.colorScheme === 'rainbow') {
      // Rainbow scheme
      if (conservation < 0.2) {
        barColor = '#E74C3C'; // Red
      } else if (conservation < 0.4) {
        barColor = '#FF7F00'; // Orange
      } else if (conservation < 0.6) {
        barColor = '#FFFF00'; // Yellow
      } else if (conservation < 0.8) {
        barColor = '#00FF00'; // Green
      } else {
        barColor = '#0000FF'; // Blue
      }
    } else if (this.colorScheme === 'heatmap') {
      // Heatmap scheme - shades of red to white
      const intensity = Math.round(conservation * 255);
      barColor = `rgb(255, ${intensity}, ${intensity})`;
    }

    const barHeight = conservation * barplotHeight;
    this.ctx.fillStyle = barColor;
    this.ctx.fillRect(
      x + barplotPadding,
      barplotTop + barplotHeight - barHeight,
      cellWidth - barplotPadding * 2,
      barHeight
    );

    // Add outline to the bar
    this.ctx.strokeStyle = 'rgba(100, 100, 100, 0.3)';
    this.ctx.lineWidth = 1;
    this.ctx.strokeRect(
      x + barplotPadding,
      barplotTop,
      cellWidth - barplotPadding * 2,
      barplotHeight
    );

    // Add conservation value text if cell is wide enough
    if (cellWidth > 20) {
      this.ctx.fillStyle = '#333333';
      this.ctx.font = '9px monospace';
      this.ctx.textAlign = 'center';
      this.ctx.textBaseline = 'middle';
      const percentText = Math.round(conservation * 100) + '%';
      this.ctx.fillText(
        percentText,
        cellCenterX,
        barplotTop + barplotHeight / 2
      );
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

  // Fixed layout properties
  private dottedCellHeight = 30; // Fixed height for dotted cells at the bottom
  private sliderHeight = 8; // Height of the slider
  private trackGap = 4; // Gap between tracks

  // Tracks
  private tracks: Map<string, MSAHeaderTrack> = new Map();

  /**
   * Constructor for the MSA Header
   * @param {MSAHeaderOptions} options - Configuration options
   */
  constructor(options: MSAHeaderOptions) {
    // Default configuration with required fields - No track-specific properties
    this.config = {
      x: options.x || 0,
      y: options.y || 0,
      width: options.width || 0,
      height: options.height || 0,
      windowStartPosition: options.windowStartPosition || 1,
      positionWidth: options.positionWidth || 15,
      totalPositions: options.totalPositions || 5000,
      headerHeight: options.headerHeight || 50,
      sliderHeight: options.sliderHeight || 8,
      currentPosition: options.currentPosition || 1,
      cellBackground: options.cellBackground !== undefined ? options.cellBackground : true,
      sliderColor: options.sliderColor || 'rgba(220, 220, 220, 0.4)',
      onPositionChange: options.onPositionChange || ((_, __) => { }),
      ...options // Override defaults with any provided options
    };

    this.sliderHeight = this.config.sliderHeight;

    // Create event element
    this.eventElement = ui.div();
    this.eventElement.style.position = 'absolute';
    this.config.canvas.parentElement?.appendChild(this.eventElement);

    // Internal state
    this.state = {
      isDragging: false,
      dragStartX: 0
    };

    // Set up event listeners
    this.setupEventListeners();
    this.init();
  }

  /**
   * Set up event listeners
   */
  private setupEventListeners(): void {
    this.eventElement.addEventListener('mousemove', (e) => {
      if (!this.isValid) return;
      if (this.isInSliderDraggableArea(e)) {
        this.eventElement.style.cursor = 'grab';
      } else if (this.isInSliderArea(e)) {
        this.eventElement.style.cursor = 'pointer';
      } else if (this.isInHeaderArea(e)) {
        this.eventElement.style.cursor = 'pointer';
      } else {
        this.eventElement.style.cursor = 'default';
      }
    });

    this.eventElement.addEventListener('mousedown', this.handleMouseDown.bind(this));
    this.eventElement.addEventListener('mousemove', this.handleMouseMove.bind(this));
    this.eventElement.addEventListener('mouseup', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('mouseleave', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('click', this.handleClick.bind(this));
    this.eventElement.addEventListener('wheel', this.handleMouseWheel.bind(this));
    window.addEventListener('keydown', this.handleKeyDown.bind(this));
  }

  /**
   * Initialize the component
   */
  private init(): void {
    // Get canvas and context
    this.canvas = this.config.canvas;
    if (!this.canvas) {
      console.error(`canvas not found.`);
      return;
    }

    const context = this.canvas.getContext('2d');
    if (!context) {
      console.error('Failed to get 2D context from canvas');
      return;
    }
    this.ctx = context;

    // Initialize all tracks
    this.tracks.forEach(track => track.init(context));
  }

  /**
   * Add a new track to the header
   */
  public addTrack(id: string, track: MSAHeaderTrack): void {
    if (this.ctx) {
      track.init(this.ctx);
    }
    this.tracks.set(id, track);
  }

  /**
   * Remove a track from the header
   */
  public removeTrack(id: string): void {
    this.tracks.delete(id);
  }

  /**
   * Get a track by ID
   */
  public getTrack<T extends MSAHeaderTrack>(id: string): T | undefined {
    return this.tracks.get(id) as T | undefined;
  }

  /**
   * Update a track's data
   * @param id The track ID
   * @param updater A function that receives the track and updates it
   */
  public updateTrack<T extends MSAHeaderTrack>(id: string, updater: (track: T) => void): void {
    const track = this.getTrack<T>(id);
    if (track) updater(track);
  }

  public get isValid() {
    return !!this.canvas && !!this.ctx && this.config.height >= this.dottedCellHeight + this.sliderHeight;
  }

  /**
   * Calculate total height needed for all visible tracks
   */
  private calculateTotalTracksHeight(): number {
    let totalHeight = 0;
    let firstTrack = true;

    // Calculate height for all visible tracks
    this.tracks.forEach(track => {
      if (track.isVisible()) {
        if (!firstTrack) {
          totalHeight += this.trackGap;
        } else {
          firstTrack = false;
        }
        totalHeight += track.getHeight();
      }
    });

    return totalHeight;
  }

  /**
   * Determine visible tracks based on available height
   */
  private determineVisibleTracks(): void {
    // Get available height for tracks (excluding dotted cells and slider and a gap)
    const availableHeight = this.config.headerHeight - (this.dottedCellHeight + this.sliderHeight + this.trackGap);

    if (availableHeight <= 0) {
      // Hide all tracks if there's not enough space
      this.tracks.forEach(track => track.setVisible(false));
      return;
    }

    // Calculate how much space we need for all tracks including gaps
    let requiredHeight = 0;
    let trackCount = 0;

    // First pass: calculate total height needed
    this.tracks.forEach(track => {
      if (trackCount > 0) requiredHeight += this.trackGap; // Add gap for each track after the first
      requiredHeight += track.getMinHeight();
      trackCount++;
    });

    if (requiredHeight <= availableHeight) {
      // We have enough space for all tracks
      this.tracks.forEach(track => track.setVisible(true));
    } else {
      // We don't have enough space, need to prioritize
      // Use the order of tracks in the map for priority
      // (Last added tracks get higher priority)

      // First pass: hide all tracks
      this.tracks.forEach(track => track.setVisible(false));

      // Second pass: show tracks in reverse order until we run out of space
      // This ensures that tracks added later (higher priority in our implementation)
      // remain visible longer when space is constrained
      let remainingHeight = availableHeight;
      let isFirstTrack = true;

      // Convert to array, reverse it (to start with highest priority), and then iterate
      const tracksArray = Array.from(this.tracks.entries()).reverse();

      for (const [id, track] of tracksArray) {
        const trackHeight = track.getMinHeight();
        const spaceNeeded = isFirstTrack ? trackHeight : trackHeight + this.trackGap;

        if (remainingHeight >= spaceNeeded) {
          track.setVisible(true);
          remainingHeight -= spaceNeeded;
          isFirstTrack = false;
        } else {
          break;
        }
      }
    }
  }

  public draw(x: number, y: number, w: number, h: number, currentPos: number, scrollerStart: number, preventable: Preventable): void {
    // Update internal state
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

    this.ctx!.save();
    // Clear canvas
    this.ctx!.clearRect(x, y, w, h);
    this.ctx!.translate(x, y);
    this.ctx!.rect(0, 0, w, h);
    this.ctx!.clip();

    // Determine which tracks should be visible based on available height
    this.determineVisibleTracks();

    // Calculate positions for tracks and dotted cells
    const sliderTop = h - this.sliderHeight;
    const dottedCellsTop = sliderTop - this.dottedCellHeight;

    // Calculate total height needed for tracks
    const totalTracksHeight = this.calculateTotalTracksHeight();

    // Position for first track
    let trackY = Math.max(0, dottedCellsTop - totalTracksHeight - this.trackGap); // Add gap between tracks and dotted cells

    // Collect visible tracks to draw connections
    const visibleTrackPositions: { y: number, height: number }[] = [];

    // Draw all visible tracks
    this.tracks.forEach(track => {
      if (track.isVisible()) {
        const trackHeight = track.getHeight();

        track.draw(
          0,
          trackY,
          w,
          trackHeight,
          this.config.windowStartPosition,
          this.config.positionWidth,
          this.config.totalPositions,
          this.config.currentPosition
        );

        // Store track position and height for drawing connections later
        visibleTrackPositions.push({ y: trackY, height: trackHeight });

        trackY += trackHeight + this.trackGap;
      }
    });

    // Draw dotted cells (always at the bottom)
    this.drawDottedCells(0, dottedCellsTop, w, this.dottedCellHeight, sliderTop);

    // Add dotted cells to track positions for connections
    visibleTrackPositions.push({ y: dottedCellsTop, height: this.dottedCellHeight });

    // Draw connection lines between tracks for selected position
    if (this.config.currentPosition >= 1 &&
      this.config.currentPosition <= this.config.totalPositions) {

      const cellWidth = this.config.positionWidth;
      const position = this.config.currentPosition;
      const windowStart = this.config.windowStartPosition;

      // Calculate position to draw connector
      const visibleIndex = position - windowStart;
      if (visibleIndex >= 0 && visibleIndex < Math.floor(w / cellWidth)) {
        const posX = visibleIndex * cellWidth;
        const cellCenterX = posX + cellWidth / 2;

        // Draw connecting lines between all visible tracks
        for (let i = 0; i < visibleTrackPositions.length - 1; i++) {
          const upperTrack = visibleTrackPositions[i];
          const lowerTrack = visibleTrackPositions[i + 1];

          // Draw connecting line
          this.ctx!.strokeStyle = 'rgba(60, 177, 115, 0.4)'; // Dim green, same as selection
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
   * Draw the dotted cells area
   */
  private drawDottedCells(x: number, y: number, width: number, height: number, sliderTop: number): void {
    if (!this.ctx) return;

    const totalPositions = this.config.totalPositions;
    const positionWidth = this.config.positionWidth;
    const currentPosition = this.config.currentPosition;
    const windowStart = this.config.windowStartPosition;
    const visiblePositionsN = Math.floor(width / positionWidth);
    const topPadding = 5;
    const posIndexTop = y + topPadding;

    // Draw the slider
    this.drawSlider(x, sliderTop, width);

    // Draw position dots and numbers
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > totalPositions) break;

      const posX = x + (i * positionWidth);
      const cellWidth = positionWidth;
      const cellCenterX = posX + cellWidth / 2;

      // Draw cell background
      if (this.config.cellBackground) {
        // Alternating cell background
        this.ctx.fillStyle = i % 2 === 0 ? 'rgba(248, 248, 248, 0.3)' : 'rgba(242, 242, 242, 0.2)';
        this.ctx.fillRect(posX, y, cellWidth, height);

        // Cell borders
        this.ctx.strokeStyle = 'rgba(220, 220, 220, 0.7)';
        this.ctx.beginPath();
        this.ctx.moveTo(posX, y);
        this.ctx.lineTo(posX, sliderTop);
        this.ctx.stroke();
      }

      // Draw position dot
      this.ctx.fillStyle = '#999999';
      this.ctx.beginPath();
      this.ctx.arc(cellCenterX, posIndexTop + 5, 1, 0, Math.PI * 2);
      this.ctx.fill();

      // Draw position number for every 10th position or current position
      if (position === currentPosition || ((position === 1 || position % 10 === 0) && Math.abs(position - currentPosition) > 1)) {
        this.ctx.fillStyle = '#333333';
        this.ctx.font = '12px monospace';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText(position.toString(), cellCenterX, posIndexTop + 15);
      }

      // Highlight current selected position
      if (position === currentPosition) {
        this.ctx.fillStyle = 'rgba(60, 177, 115, 0.2)';
        this.ctx.fillRect(
          posX,
          y,
          cellWidth,
          height
        );
      }
    }
  }

  /**
   * Draw the slider
   */
  private drawSlider(x: number, sliderTop: number, width: number): void {
    if (!this.ctx) return;

    // Draw the full sequence slider bar
    this.ctx.fillStyle = this.config.sliderColor;
    this.ctx.fillRect(x, sliderTop, width, this.sliderHeight);

    const visiblePositionsN = Math.floor(width / this.config.positionWidth);
    const windowStart = this.config.windowStartPosition;

    // Calculate slider position on the bar
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderWidth = this.sliderWidth;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      (windowStart - 1) / totalSliderRange * (width - sliderWidth);

    const sliderLengthPX = totalSliderRange <= 0 ? width : sliderWidth;

    // Draw slider window
    this.ctx.fillStyle = 'rgba(150, 150, 150, 0.5)';
    this.ctx.fillRect(x + sliderStartPX, sliderTop, sliderLengthPX, this.sliderHeight);

    // Draw position marker if current position is valid
    if (this.config.currentPosition >= 1 && this.config.currentPosition <= this.config.totalPositions) {
      // Calculate position on the bar
      const currentPositionRatio = (this.config.currentPosition - 1) / (this.config.totalPositions - 1);
      const notchX = Math.round(currentPositionRatio * width);

      // Draw marker
      this.ctx.fillStyle = '#3CB173'; // Green
      this.ctx.fillRect(
        x + notchX - 1,
        sliderTop - 2,
        3,
        this.sliderHeight + 4
      );
    }
  }

  /**
   * Set up the event element
   */
  private setupEventElement(): void {
    this.eventElement.style.display = 'block';
    this.eventElement.style.left = `${this.config.x}px`;
    this.eventElement.style.top = `${this.config.y}px`;
    this.eventElement.style.width = `${this.config.width}px`;
    this.eventElement.style.height = `${this.config.height}px`;
  }

  /**
   * Get coordinates relative to the header
   */
  getCoords(e: MouseEvent) {
    const rect = this.canvas!.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;
    const y = e.clientY - rect.top - this.config.y;
    return { x, y };
  }

  isInHeaderArea(e: MouseEvent): boolean {
    const { x, y } = this.getCoords(e);
    return x >= 0 && x <= this.config.width && y >= 0 && y <= this.config.headerHeight;
  }

  get positionWidth(): number {
    return this.config.positionWidth;
  }

  public set positionWidth(value: number) {
    this.config.positionWidth = value;
  }

  isInSliderArea(e: MouseEvent): boolean {
    const { y } = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.sliderHeight;
    return y > sliderTop && y < sliderTop + this.sliderHeight;
  }

  get sliderWidth(): number {
    const pseudoPositionWidth = this.config.width / this.config.totalPositions;
    const w = pseudoPositionWidth * (this.config.width / this.config.positionWidth);
    return Math.max(w, 20);
  }

  isInSliderDraggableArea(e: MouseEvent): boolean {
    const { x, y } = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.sliderHeight;
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = this.config.windowStartPosition;

    // Calculate slider position
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      (windowStart - 1) / totalSliderRange * (this.config.width - this.sliderWidth);

    return y > sliderTop && y < sliderTop + this.sliderHeight &&
      x >= sliderStartPX && x < sliderStartPX + this.sliderWidth;
  }

  /**
   * Handle mouse down (start dragging)
   */
  private handleMouseDown(e: MouseEvent): void {
    if (!this.isValid) return;
    const { x } = this.getCoords(e);
    if (this.isInSliderDraggableArea(e)) {
      this.state.isDragging = true;
      this.state.dragStartX = x;
      this.handleSliderDrag(x);
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
    }
  }

  /**
   * Handle mouse wheel
   */
  private handleMouseWheel(e: WheelEvent): void {
    if (!this.isValid) return;
    if (this.isInHeaderArea(e)) {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      // Determine scroll direction and amount
      const delta = e.shiftKey ? Math.sign(e.deltaY) : Math.sign(e.deltaX || e.deltaY);
      const scrollSpeed = e.shiftKey ? 3 : 1; // Faster scrolling with shift
      const newStartPosition = this.config.windowStartPosition + (delta * scrollSpeed);

      // Clamp to valid range
      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const maxStart = this.config.totalPositions - visiblePositions + 1;
      this.config.windowStartPosition = Math.max(1, Math.min(maxStart, newStartPosition));

      if (typeof this.config.onPositionChange === 'function') {
        this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
      }
    }
  }

  /**
   * Handle mouse move (dragging)
   */
  private handleMouseMove(e: MouseEvent): void {
    if (!this.state.isDragging || !this.isValid) return;

    const rect = this.canvas!.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;

    this.handleSliderDrag(x);
    e.preventDefault();
    e.stopPropagation();
    e.stopImmediatePropagation();
  }

  /**
   * Handle keyboard navigation
   */
  private handleKeyDown(e: KeyboardEvent): void {
    if (!this.isValid || this.config.currentPosition < 1) return;
    if (!document.activeElement?.contains(this.eventElement) ||
      this.eventElement.style.display !== 'block') return;

    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      // Determine scroll direction
      const delta = e.key === 'ArrowLeft' ? -1 : 1;
      const newPosition = Math.min(Math.max(this.config.currentPosition + delta, 1), this.config.totalPositions);

      if (newPosition === this.config.currentPosition) return;
      this.config.currentPosition = newPosition;

      // Make sure the current position is visible
      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const start = this.config.windowStartPosition;
      const end = start + visiblePositions - 1;

      if (newPosition < start || newPosition > end) {
        if (delta < 0) {
          this.config.windowStartPosition = newPosition;
        } else {
          this.config.windowStartPosition = Math.max(1, newPosition - visiblePositions + 1);
        }
      }
    } else if (e.key === 'Escape') {
      // Reset the current position
      this.config.currentPosition = -2;
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
    } else {
      return;
    }

    if (typeof this.config.onPositionChange === 'function') {
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  /**
   * Handle mouse up (end dragging)
   */
  private handleMouseUp(): void {
    this.state.isDragging = false;
  }

  /**
   * Handle slider drag
   */
  private handleSliderDrag(x: number): void {
    if (!this.isValid) return;

    const sliderWidth = this.sliderWidth;
    const canvasWidth = this.config.width - sliderWidth;
    const normalizedX = Math.max(0, Math.min(this.config.width, x));
    const fittedPositions = Math.floor(this.config.width / this.config.positionWidth);
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;

    // Calculate new window start position
    const sliderStartPx = Math.max(0, normalizedX - sliderWidth / 2);
    const windowStart = sliderStartPx / (canvasWidth) * (totalSliderRange);
    this.config.windowStartPosition = Math.max(1, Math.min(windowStart, this.config.totalPositions - fittedPositions + 1));

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function') {
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  get headerHeight(): number {
    return this.config.headerHeight;
  }

  /**
   * Set header height without redrawing
   */
  set headerHeight(value: number) {
    this.config.headerHeight = value;
  }

  /**
   * Handle click on positions
   */
  private handleClick(e: MouseEvent): void {
    if (!this.isValid) return;

    const { x, y } = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.sliderHeight;

    if (y < sliderTop && y >= 0) {
      // Calculate which position was clicked
      const cellWidth = this.config.positionWidth;
      const clickedCellIndex = Math.round(x / cellWidth - 0.5);
      const windowStart = this.config.windowStartPosition;
      const clickedPosition = windowStart + clickedCellIndex;

      // Update current position if valid
      if (clickedPosition >= 1 && clickedPosition <= this.config.totalPositions) {
        this.config.currentPosition = clickedPosition;

        // Call callback if defined
        if (typeof this.config.onPositionChange === 'function') {
          this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
        }
      }
    }
  }

  /**
   * Get the current window range
   */
  public getWindowRange(): WindowRange {
    return {
      start: this.config.windowStartPosition,
      end: Math.min(
        this.config.totalPositions,
        this.config.windowStartPosition + Math.floor(this.config.width / this.config.positionWidth)
      )
    };
  }

  /**
   * Update configuration
   */
  public updateConfig(newConfig: Partial<MSAHeaderOptions>): void {
    // Update config with new values
    Object.assign(this.config, newConfig);

    // Ensure current position is still valid
    this.config.currentPosition = Math.min(this.config.currentPosition, this.config.totalPositions);

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function') {
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  /**
   * Get current position
   */
  public getCurrentPosition(): number {
    return this.config.currentPosition;
  }

  /**
   * Set current position
   */
  public setCurrentPosition(position: number): void {
    // Clamp to valid range
    this.config.currentPosition = Math.max(1, Math.min(this.config.totalPositions, position));

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function') {
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }
}
