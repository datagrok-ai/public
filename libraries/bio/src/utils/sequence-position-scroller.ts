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

export class MSAScrollingHeader {
  private config: Required<MSAHeaderOptions>;
  private state: MSAHeaderState;
  private canvas: HTMLCanvasElement | null = null;
  private ctx: CanvasRenderingContext2D | null = null;
  private eventElement: HTMLDivElement;

  private titlePadding = 8; // Space between title and dot cells
  private maxDotCellHeight = 30; // Fixed height for dotted cells
  private minTitleHeight = 20; // Minimum height for title
  private minDotCellsHeight = 40; // Minimum height for dotted cells section (includes slider)
  private minConservationHeight = 40; // Minimum height for conservation barplots

  // Thresholds for different rendering modes
  private titleOnlyThreshold = this.minTitleHeight;
  private dotCellsThreshold = this.minTitleHeight + this.minDotCellsHeight;
  private fullDisplayThreshold = this.minTitleHeight + this.minDotCellsHeight + this.minConservationHeight;

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
      // Add default values for conservation properties
      conservationData: options.conservationData || [],
      conservationHeight: options.conservationHeight || 40,
      conservationColorScheme: options.conservationColorScheme || 'default',
      ...options // Override defaults with any provided options
    };
    this.eventElement = ui.div();
    this.eventElement.style.position = 'absolute';
    this.config.canvas.parentElement?.appendChild(this.eventElement);

    // Internal state
    this.state = {
      isDragging: false,
      dragStartX: 0
    };

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

    this.init();
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

    // Add event listeners
    this.eventElement.addEventListener('mousedown', this.handleMouseDown.bind(this));
    this.eventElement.addEventListener('mousemove', this.handleMouseMove.bind(this));
    this.eventElement.addEventListener('mouseup', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('mouseleave', this.handleMouseUp.bind(this));
    this.eventElement.addEventListener('click', this.handleClick.bind(this));
    this.eventElement.addEventListener('wheel', this.handleMouseWheel.bind(this));
    window.addEventListener('keydown', this.handleKeyDown.bind(this));
  }

  public get isValid() {
    return !!this.canvas && !!this.ctx && this.config.height >= this.config.headerHeight;
  }

  public draw(x: number, y: number, w: number, h: number, currentPos: number, scrollerStart: number, preventable: Preventable): void {
    // soft internal update
    if (this.config.x != x || this.config.y != y || this.config.width != w || this.config.height != h || this.config.currentPosition != currentPos || this.config.windowStartPosition != scrollerStart)
      Object.assign(this.config, { x, y, width: w, height: h, currentPosition: currentPos, windowStartPosition: scrollerStart });

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

    // Calculate dimensions
    const canvasWidth = w;
    const canvasHeight = h;

    // Determine the rendering mode based on available height
    const availableHeight = this.config.headerHeight;

    // Title-only mode
    if (availableHeight <= this.titleOnlyThreshold) {
      this.drawTitle(canvasWidth);
      this.ctx!.restore();
      preventable.preventDefault();
      this.setupEventElement();
      return;
    }

    // Dotted cells mode (no conservation)
    const showConservation = availableHeight > this.dotCellsThreshold + this.minConservationHeight &&
      this.config.conservationData.length > 0 &&
      this.config.conservationHeight > 0;

    // Calculate layout based on available components
    let titleHeight = 0;
    const dottedCellsHeight = this.maxDotCellHeight;
    const conservationHeight = showConservation ? this.config.conservationHeight : 0;
    let remainingHeight = availableHeight;

    // If there's enough space for a title, allocate space for it
    if (availableHeight > this.dotCellsThreshold + this.minTitleHeight) {
      titleHeight = this.minTitleHeight;
      remainingHeight -= titleHeight;
    }

    // Calculate slider position
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;

    // Calculate dotted cells position
    const dottedCellsTop = sliderTop - dottedCellsHeight;
    const topPadding = 5;
    const posIndexTop = dottedCellsTop + topPadding;

    // Draw title if there's space
    if (titleHeight > 0) {
      this.drawTitle(canvasWidth);
    }

    // Conservation barplot dimensions
    const barplotHeight = conservationHeight;
    const barplotTop = dottedCellsTop - barplotHeight - 5; // Position barplots above dotted cells with small gap
    const barplotPadding = 2; // Padding between bars

    // Draw the full sequence slider bar
    this.ctx!.fillStyle = this.config.sliderColor;
    this.ctx!.fillRect(0, sliderTop, canvasWidth, this.config.sliderHeight);

    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = Math.max(1, this.config.windowStartPosition);

    // Calculate slider position on the bar
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      windowStart / (totalSliderRange) * (canvasWidth - this.sliderWidth);

    const sliderLengthPX = totalSliderRange <= 0 ? canvasWidth :
      this.sliderWidth;

    // Draw slider window (darker rectangle)
    this.ctx!.fillStyle = 'rgba(150, 150, 150, 0.5)';
    this.ctx!.fillRect(sliderStartPX, sliderTop, sliderLengthPX, this.config.sliderHeight);

    if (this.config.currentPosition >= 1 && this.config.currentPosition <= this.config.totalPositions) {
      // Calculate the position of the current selection as a proportion of total sequence
      const currentPositionRatio = (this.config.currentPosition - 1) / (this.config.totalPositions - 1);
      const notchX = Math.round(currentPositionRatio * canvasWidth);

      // Draw a 3px thick vertical green marker at this position
      this.ctx!.fillStyle = '#3CB173'; // Same green as the cell highlight
      this.ctx!.fillRect(
        notchX - 1, // Center the 3px marker on the calculated position
        sliderTop - 2, // Extend the marker 2px above the scrollbar
        3, // 3px width as requested
        this.config.sliderHeight + 4 // Make the marker taller than the scrollbar
      );
    }

    // Draw position marks and indices
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > this.config.totalPositions) break;

      const x = i * this.config.positionWidth;
      const cellWidth = this.config.positionWidth;
      const cellCenterX = x + cellWidth / 2;

      // Draw conservation barplot if enabled
      if (showConservation && position - 1 < this.config.conservationData.length) {
        this.drawConservationBar(
          position - 1,
          x,
          cellWidth,
          cellCenterX,
          barplotTop,
          barplotHeight,
          barplotPadding
        );
      }

      // Draw cell background for monospace appearance
      if (this.config.cellBackground) {
        // Very light alternating cell background
        this.ctx!.fillStyle = i % 2 === 0 ? 'rgba(248, 248, 248, 0.3)' : 'rgba(242, 242, 242, 0.2)';
        this.ctx!.fillRect(x, dottedCellsTop, cellWidth, dottedCellsHeight);

        // Cell borders - very light vertical lines
        this.ctx!.strokeStyle = 'rgba(220, 220, 220, 0.7)';
        this.ctx!.beginPath();
        this.ctx!.moveTo(x, dottedCellsTop);
        this.ctx!.lineTo(x, sliderTop);
        this.ctx!.stroke();
      }

      // Draw position dot for every position - centered in cell
      this.ctx!.fillStyle = '#999999';
      this.ctx!.beginPath();
      this.ctx!.arc(cellCenterX, posIndexTop + 5, 1, 0, Math.PI * 2);
      this.ctx!.fill();

      // Draw position number for every 10th position, on the current position and also make sure that the number is not obstructed by some other numbers
      if (position === this.config.currentPosition || ((position === 1 || position % 10 === 0) && Math.abs(position - this.config.currentPosition) > 1)) {
        this.ctx!.fillStyle = '#333333';
        this.ctx!.font = '12px monospace';
        this.ctx!.textAlign = 'center';
        this.ctx!.textBaseline = 'middle';
        this.ctx!.fillText(position.toString(), cellCenterX, posIndexTop + 15);
      }

      // Highlight current selected position with square marker
      if (position === this.config.currentPosition) {
        // Draw filled background with 50% opacity
        this.ctx!.fillStyle = 'rgba(60, 177, 115, 0.2)';
        this.ctx!.fillRect(
          x,
          dottedCellsTop,
          cellWidth,
          dottedCellsHeight
        );

        // Also highlight the conservation bar for the current position
        if (showConservation && position - 1 < this.config.conservationData.length) {
          this.ctx!.fillStyle = 'rgba(60, 177, 115, 0.1)';
          this.ctx!.fillRect(
            x,
            barplotTop,
            cellWidth,
            barplotHeight
          );

          // Add a vertical line connecting the highlighted areas
          this.ctx!.strokeStyle = 'rgba(60, 177, 115, 0.4)';
          this.ctx!.lineWidth = 1;
          this.ctx!.beginPath();
          this.ctx!.moveTo(x + cellWidth / 2, barplotTop + barplotHeight);
          this.ctx!.lineTo(x + cellWidth / 2, dottedCellsTop);
          this.ctx!.stroke();
        }
      }
    }

    this.ctx!.restore();
    // Prevent default behavior if in header area
    preventable.preventDefault();
    this.setupEventElement();
  }
  private setupEventElement(): void {
    this.eventElement.style.display = 'block';
    this.eventElement.style.left = `${this.config.x}px`;
    this.eventElement.style.top = `${this.config.y}px`;
    this.eventElement.style.width = `${this.config.width}px`;
    this.eventElement.style.height = `${this.config.height}px`;
  }
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

    const conservation = this.config.conservationData[posIndex];

    // Draw bar background
    this.ctx.fillStyle = 'rgba(240, 240, 240, 0.5)';
    this.ctx.fillRect(
      x + barplotPadding,
      barplotTop,
      cellWidth - barplotPadding * 2,
      barplotHeight
    );

    // Draw conservation bar with color based on color scheme
    let barColor = '#3CB173'; // Default green for high conservation

    if (this.config.conservationColorScheme === 'default') {
      // Default scheme: green (high), yellow (medium), red (low)
      if (conservation < 0.5) {
        barColor = '#E74C3C'; // Red for low conservation (<50%)
      } else if (conservation < 0.75) {
        barColor = '#F39C12'; // Yellow for medium conservation (50-75%)
      }
    } else if (this.config.conservationColorScheme === 'rainbow') {
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
    } else if (this.config.conservationColorScheme === 'heatmap') {
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

  private drawTitle(width: number): void {
    if (!this.ctx) return;

    this.ctx.font = 'bold 13px Roboto, Roboto Local';
    this.ctx.fillStyle = '#4a4a49';
    this.ctx.textAlign = 'center';
    this.ctx.textBaseline = 'middle';

    // Position the title in the middle of the available space
    const titleX = width / 2;
    const titleY = this.minTitleHeight / 2;

    // Draw the text - if we're in a grid context, this would be the column name
    this.ctx.fillText('Sequence', titleX, titleY);
  }
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
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;
    return y > sliderTop && y < sliderTop + this.config.sliderHeight;
  }

  get sliderWidth(): number {
    const pseudoPositionWidth = this.config.width / this.config.totalPositions;
    const w = pseudoPositionWidth * (this.config.width / this.config.positionWidth);
    return Math.max(w, 20);
  }

  isInSliderDraggableArea(e: MouseEvent): boolean {
    const { x, y } = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = Math.max(1, this.config.windowStartPosition);
    // Calculate slider position on the bar
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      windowStart / (totalSliderRange) * (this.config.width - this.sliderWidth);


    return y > sliderTop && y < sliderTop + this.config.sliderHeight && x >= sliderStartPX && x < sliderStartPX + this.sliderWidth;
  }

  /**
   * Handle mouse down (start dragging)
   * @param {MouseEvent} e - Mouse event
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

  private handleMouseWheel(e: WheelEvent): void {
    if (!this.isValid) return;
    if (this.isInHeaderArea(e)) {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      // Determine scroll direction and amount
      // Use deltaX if available (for horizontal scrolling devices)
      // otherwise use deltaY with shift key for horizontal scroll
      const delta = e.shiftKey ? Math.sign(e.deltaY) : Math.sign(e.deltaX || e.deltaY);

      // Move by 1 position or more for faster scrolling (optional)
      const scrollSpeed = e.shiftKey ? 3 : 1; // Faster scrolling with shift
      const newStartPosition = this.config.windowStartPosition + (delta * scrollSpeed);

      // Clamp to valid range
      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const maxStart = this.config.totalPositions - visiblePositions + 1;
      this.config.windowStartPosition = Math.max(1, Math.min(maxStart, newStartPosition));

      if (typeof this.config.onPositionChange === 'function')
        this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  /**
   * Handle mouse move (dragging)
   * @param {MouseEvent} e - Mouse event
   */
  private handleMouseMove(e: MouseEvent): void {
    if (!this.state.isDragging || !this.isValid)
      return;

    const rect = this.canvas!.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;

    this.handleSliderDrag(x);
    e.preventDefault();
    e.stopPropagation();
    e.stopImmediatePropagation();
  }

  private handleKeyDown(e: KeyboardEvent): void {
    if (!this.isValid || this.config.currentPosition < 1) return;
    if (!document.activeElement?.contains(this.eventElement) || this.eventElement.style.display !== 'block')
      return;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();

      // Determine scroll direction
      const delta = e.key === 'ArrowLeft' ? -1 : 1;

      const newPosition = Math.min(Math.max(this.config.currentPosition + delta, 1), this.config.totalPositions);
      if (newPosition === this.config.currentPosition)
        return;
      this.config.currentPosition = newPosition;
      // make sure that the cureent position is visible
      const visiblePositions = Math.floor(this.config.width / this.config.positionWidth);
      const start = this.config.windowStartPosition;
      const end = start + visiblePositions - 1;
      if (newPosition < start || newPosition > end) {
        if (delta < 0)
          this.config.windowStartPosition = newPosition;
        else
          this.config.windowStartPosition = Math.max(1, newPosition - visiblePositions + 1);
      }
    } else if (e.key === 'Escape') {
      // reset the current position
      this.config.currentPosition = -2;
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
    } else return;

    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  /**
   * Handle mouse up (end dragging)
   */
  private handleMouseUp(): void {
    this.state.isDragging = false;
  }

  /**
   * Handle slider drag
   * @param {number} x - X position of mouse
   */
  private handleSliderDrag(x: number): void {
    if (!this.isValid) return;

    const sliderWidth = this.sliderWidth;
    // Calculate the position based on the drag position
    const canvasWidth = this.config.width - sliderWidth;

    const normalizedX = Math.max(0, Math.min(this.config.width, x));
    const fittedPositions = Math.floor(this.config.width / this.config.positionWidth);

    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    // Calculate slider position on the bar
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    // const sliderStartPX =
    //   windowStart / (totalSliderRange) * (this.config.width - this.sliderWidth);
    // after we normalize the x position of the mouse, this is where the center of the slider should be.
    const sliderStartPx = Math.max(0, normalizedX - sliderWidth / 2);
    // then we reverse the formula and calculate the new start position
    const windowStart = sliderStartPx / (canvasWidth) * (totalSliderRange);
    // we add these positions so that it feels like we are grabbing the slider by its center
    this.config.windowStartPosition = Math.max(1, Math.min(windowStart, this.config.totalPositions - fittedPositions + 1));

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  get headerHeight(): number {
    return this.config.headerHeight;
  }

  /** Soft setting of header hight, without redrawing or update event
   * @param {number} value - New header height
   */
  set headerHeight(value: number) {
    this.config.headerHeight = value;
  }

  /**
   * Handle click on positions
   * @param {MouseEvent} e - Mouse event
   */
  private handleClick(e: MouseEvent): void {
    if (!this.isValid) return;

    // Get calculated coordinates
    const { x, y } = this.getCoords(e);


    const sliderTop = this.config.headerHeight - this.config.sliderHeight;
    if (y < sliderTop && y >= 0) {
      // Calculate which position was clicked
      const cellWidth = this.config.positionWidth;
      const clickedCellIndex = Math.round(x / cellWidth - 0.5);

      // Calculate actual position in sequence
      const windowStart = this.config.windowStartPosition;
      const clickedPosition = windowStart + clickedCellIndex;

      // Update current position if valid
      if (clickedPosition >= 1 && clickedPosition <= this.config.totalPositions) {
        this.config.currentPosition = clickedPosition;

        // Call callback if defined
        if (typeof this.config.onPositionChange === 'function')
          this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
      }
    }
  }

  /**
   * Get the current window range
   * @return {WindowRange} Object with start and end properties
   */
  public getWindowRange(): WindowRange {
    return {
      start: this.config.windowStartPosition,
      end: Math.min(this.config.totalPositions, this.config.windowStartPosition + Math.floor(this.config.width / this.config.positionWidth))
    };
  }

  /**
   * Update configuration
   * @param {MSAHeaderOptions} newConfig - New configuration options
   */
  public updateConfig(newConfig: Partial<MSAHeaderOptions>): void {
    // Update config with new values
    Object.assign(this.config, newConfig);

    // Ensure current position is still valid
    this.config.currentPosition = Math.min(this.config.currentPosition, this.config.totalPositions);

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }

  /**
   * Get current position
   * @return {number} Current position
   */
  public getCurrentPosition(): number {
    return this.config.currentPosition;
  }

  /**
   * Set current position
   * @param {number} position - New position
   */
  public setCurrentPosition(position: number): void {
    // Clamp to valid range
    this.config.currentPosition = Math.max(1, Math.min(this.config.totalPositions, position));

    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
  }
}
