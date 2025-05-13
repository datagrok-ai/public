/* eslint-disable max-len */
/**
 * MSAHeader.ts - Interactive MSA Header Component
 * This module provides functionality for an interactive MSA (Multiple Sequence Alignment) header
 * with position markers, slider navigation, and position selection.
 */

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

  /**
   * Constructor for the MSA Header
   * @param {MSAHeaderOptions} options - Configuration options
   */
  constructor(options: MSAHeaderOptions) {
    // Default configuration with required fields
    this.config = {
      x: options.x || 0,
      y: options.y || 0,
      width: options.width || 400,
      height: options.height || 60,
      windowStartPosition: options.windowStartPosition || 1,
      positionWidth: options.positionWidth || 15,
      totalPositions: options.totalPositions || 5000,
      headerHeight: options.headerHeight || 50,
      sliderHeight: options.sliderHeight || 8,
      currentPosition: options.currentPosition || 1,
      cellBackground: options.cellBackground !== undefined ? options.cellBackground : true,
      sliderColor: options.sliderColor || 'rgba(220, 220, 220, 0.4)',
      onPositionChange: options.onPositionChange || ((_, __) => {}),
      ...options // Override defaults with any provided options
    };

    // Internal state
    this.state = {
      isDragging: false,
      dragStartX: 0
    };

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
    this.canvas.addEventListener('mousedown', this.handleMouseDown.bind(this));
    this.canvas.addEventListener('mousemove', this.handleMouseMove.bind(this));
    this.canvas.addEventListener('mouseup', this.handleMouseUp.bind(this));
    this.canvas.addEventListener('mouseleave', this.handleMouseUp.bind(this));
    this.canvas.addEventListener('click', this.handleClick.bind(this));
  }

  public draw(x: number, y: number, w: number, h: number, currentPos: number, scrollerStart: number): void {
    if (!this.ctx || !this.canvas) return;

    if (this.config.currentPosition != currentPos || this.config.windowStartPosition != scrollerStart || this.config.x != x || this.config.y != y || this.config.width != w || this.config.height != h) {
      this.updateConfig({x, y, width: w, height: h, currentPosition: currentPos, windowStartPosition: scrollerStart});
      return;
    }

    this.ctx.save();
    // Clear canvas
    this.ctx.clearRect(x, y, w, h);
    this.ctx.translate(x, y);
    this.ctx.rect(0, 0, w, h);
    this.ctx.clip();

    // Calculate dimensions
    const canvasWidth = w;
    const canvasHeight = h;
    const topPadding = 5;
    const posIndexTop = topPadding;
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;

    // Draw the full sequence slider bar (very subtle gray bar)
    this.ctx.fillStyle = this.config.sliderColor;
    this.ctx.fillRect(0, sliderTop, canvasWidth, this.config.sliderHeight);

    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = Math.max(1, this.config.windowStartPosition);

    // Calculate slider position on the bar
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;
    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      windowStart / (this.config.totalPositions) * canvasWidth;

    const sliderLengthPX = totalSliderRange <= 0 ? canvasWidth :
      this.sliderWidth;

    // Draw slider window (darker rectangle)
    this.ctx.fillStyle = 'rgba(150, 150, 150, 0.5)';
    this.ctx.fillRect(sliderStartPX, sliderTop, sliderLengthPX, this.config.sliderHeight);

    // Draw position marks and indices
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > this.config.totalPositions) break;

      const x = i * this.config.positionWidth;
      const cellWidth = this.config.positionWidth;
      const cellCenterX = x + cellWidth / 2;

      // Draw cell background for monospace appearance
      if (this.config.cellBackground) {
        // Very light alternating cell background
        this.ctx.fillStyle = i % 2 === 0 ? 'rgba(248, 248, 248, 0.3)' : 'rgba(242, 242, 242, 0.2)';
        this.ctx.fillRect(x, posIndexTop, cellWidth, Math.min(this.config.headerHeight - topPadding, canvasHeight - posIndexTop) - this.config.sliderHeight);

        // Cell borders - very light vertical lines
        this.ctx.strokeStyle = 'rgba(220, 220, 220, 0.7)';
        this.ctx.beginPath();
        this.ctx.moveTo(x, posIndexTop);
        this.ctx.lineTo(x, Math.min(this.config.headerHeight, canvasHeight) - this.config.sliderHeight);
        this.ctx.stroke();
      }

      // Draw position dot for every position - centered in cell
      this.ctx.fillStyle = '#999999';
      this.ctx.beginPath();
      this.ctx.arc(cellCenterX, posIndexTop * 2, 1, 0, Math.PI * 2);
      this.ctx.fill();

      // Draw position number for every 10th position
      if (position % 10 === 0) {
        this.ctx.fillStyle = '#333333';
        this.ctx.font = '10px monospace';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText(position.toString(), cellCenterX, posIndexTop * 4);
      }

      // Highlight current selected position with square marker
      // if current is outside the range, draw it at the end
      if (position === this.config.currentPosition) {
        // Draw square indicator
        this.ctx.fillStyle = '#50A9C5';
        this.ctx.fillRect(
          x,
          0,
          cellWidth,
          sliderTop
        );
      }
    }

    this.ctx.restore();
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

  isInSliderArea(e: MouseEvent): boolean {
    const {x, y} = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;
    return y > sliderTop && y < sliderTop + this.config.sliderHeight;
  }

  get sliderWidth(): number {
    const pseudoPositionWidth = this.config.width / this.config.totalPositions;
    const w = pseudoPositionWidth * (this.config.width / this.config.positionWidth);
    return Math.max(w, 15);
  }

  isInSliderDraggableArea(e: MouseEvent): boolean {
    const {x, y} = this.getCoords(e);
    const sliderTop = this.config.headerHeight - this.config.sliderHeight;
    const startSeqX = this.config.windowStartPosition;
    const pseudoPositionWidth = this.config.width / this.config.totalPositions;
    const startSliderXPX = (startSeqX - 1) * pseudoPositionWidth;
    const endSliderXPX = startSliderXPX + this.sliderWidth;
    return y > sliderTop && y < sliderTop + this.config.sliderHeight && x >= startSliderXPX && x < endSliderXPX;
  }

  /**
   * Handle mouse down (start dragging)
   * @param {MouseEvent} e - Mouse event
   */
  private handleMouseDown(e: MouseEvent): void {
    if (!this.canvas) return;
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
    if (!this.canvas) return;
    if (this.isInHeaderArea(e)) {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
      const delta = Math.sign(e.deltaY);
      const newStartPosition = this.config.windowStartPosition + delta;
      this.config.windowStartPosition = Math.max(1, Math.min(this.config.totalPositions - Math.floor(this.config.width / this.config.positionWidth), newStartPosition));
      // Call callback if defined
      if (typeof this.config.onPositionChange === 'function')
        this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    }
  }

  /**
   * Handle mouse move (dragging)
   * @param {MouseEvent} e - Mouse event
   */
  private handleMouseMove(e: MouseEvent): void {
    if (!this.state.isDragging || !this.canvas)
      return;

    const rect = this.canvas.getBoundingClientRect();
    const x = e.clientX - rect.left - this.config.x;

    this.handleSliderDrag(x);
    e.preventDefault();
    e.stopPropagation();
    e.stopImmediatePropagation();
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
    if (!this.canvas) return;

    // Calculate the position based on the drag position
    const canvasWidth = this.config.width;

    const dragAmmount = x - this.state.dragStartX;
    const pseudoPositionToPixelRatio = canvasWidth / this.config.totalPositions;
    const fittedPositions = Math.ceil(canvasWidth / this.config.positionWidth);
    const newStartPosition = Math.round(this.config.windowStartPosition + dragAmmount / pseudoPositionToPixelRatio);
    this.config.windowStartPosition = Math.max(1, Math.min(this.config.totalPositions - fittedPositions, newStartPosition));
    // Call callback if defined
    if (typeof this.config.onPositionChange === 'function')
      this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    this.state.dragStartX = x;
  }

  /**
   * Handle click on positions
   * @param {MouseEvent} e - Mouse event
   */
  private handleClick(_e: MouseEvent): void {
    // if (!this.canvas) return;

    // const rect = this.canvas.getBoundingClientRect();
    // const x = e.clientX - rect.left;
    // const y = e.clientY - rect.top;

    // const sliderTop = 10;
    // const posIndexTop = sliderTop + this.config.sliderHeight + 5;

    // // Process clicks in position indicator area (not in slider area)
    // if (y > posIndexTop) {
    //   // Calculate column clicked
    //   const columnIndex = Math.floor(x / this.config.columnWidth);

    //   // Calculate the actual position in the sequence
    //   const halfVisible = Math.floor(this.config.visibleWidth / 2);
    //   const windowStart = Math.max(1, this.config.currentPosition - halfVisible);
    //   const windowEnd = Math.min(this.config.totalPositions, windowStart + this.config.visibleWidth - 1);
    //   const adjustedWindowStart = Math.max(1, windowEnd - this.config.visibleWidth + 1);

    //   const clickedPosition = adjustedWindowStart + columnIndex;

    //   // Update current position if valid
    //   if (clickedPosition >= 1 && clickedPosition <= this.config.totalPositions) {
    //     this.config.currentPosition = clickedPosition;

    //     // Call callback if defined
    //     if (typeof this.config.onPositionChange === 'function')
    //       this.config.onPositionChange(this.config.currentPosition, this.getWindowRange());
    //   }
    // }
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
