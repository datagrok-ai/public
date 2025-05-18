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
  conservationData?: number[];
  conservationHeight?: number;
  conservationColorScheme?: 'default' | 'rainbow' | 'heatmap';
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
abstract class MSAHeaderTrack {
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

/**
 * Track for displaying conservation bars
 */
class ConservationTrack extends MSAHeaderTrack {
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
  private trackGap = 5; // Gap between tracks

  // Tracks
  private tracks: Map<string, MSAHeaderTrack> = new Map();

  /**
   * Constructor for the MSA Header
   * @param {MSAHeaderOptions} options - Configuration options
   */
  constructor(options: MSAHeaderOptions) {
    // Default configuration with required fields
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
      conservationData: options.conservationData || [],
      conservationHeight: options.conservationHeight || 40,
      conservationColorScheme: options.conservationColorScheme || 'default',
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

    // Initialize tracks
    if (this.config.conservationData && this.config.conservationData.length > 0) {
      const conservationTrack = new ConservationTrack(
        this.config.conservationData,
        this.config.conservationHeight,
        this.config.conservationColorScheme
      );
      this.tracks.set('conservation', conservationTrack);
    }

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
  public getTrack(id: string): MSAHeaderTrack | undefined {
    return this.tracks.get(id);
  }

  /**
   * Update conservation data
   */
  public updateConservationData(data: number[]): void {
    let conservationTrack = this.tracks.get('conservation') as ConservationTrack;

    if (conservationTrack) {
      conservationTrack.updateData(data);
    } else if (data && data.length > 0) {
      conservationTrack = new ConservationTrack(
        data,
        this.config.conservationHeight,
        this.config.conservationColorScheme
      );
      if (this.ctx) {
        conservationTrack.init(this.ctx);
      }
      this.tracks.set('conservation', conservationTrack);
    }
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

    this.tracks.forEach(track => {
      if (trackCount > 0) requiredHeight += this.trackGap; // Add gap for each track after the first
      requiredHeight += track.getMinHeight();
      trackCount++;
    });

    if (requiredHeight <= availableHeight) {
      // We have enough space for all tracks
      this.tracks.forEach(track => track.setVisible(true));
    } else {
      // We don't have enough space, prioritize tracks
      // For now, just show the conservation track if there's enough space
      const conservationTrack = this.tracks.get('conservation');
      if (conservationTrack && availableHeight >= conservationTrack.getMinHeight()) {
        conservationTrack.setVisible(true);

        // Hide other tracks
        this.tracks.forEach((track, id) => {
          if (id !== 'conservation') {
            track.setVisible(false);
          }
        });
      } else {
        // Hide all tracks
        this.tracks.forEach(track => track.setVisible(false));
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

    // Update conservation data if provided
    if (newConfig.conservationData) {
      this.updateConservationData(newConfig.conservationData);
    }

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
