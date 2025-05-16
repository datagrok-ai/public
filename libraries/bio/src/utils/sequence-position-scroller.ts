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
      conservationHeight: options.conservationHeight || 0,
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
    // Update internal state
    if (this.config.x != x || this.config.y != y || this.config.width != w || this.config.height != h ||
      this.config.currentPosition != currentPos || this.config.windowStartPosition != scrollerStart) {
      Object.assign(this.config, { x, y, width: w, height: h, currentPosition: currentPos, windowStartPosition: scrollerStart });
    }

    if (!this.isValid) {
      this.eventElement.style.display = 'none';
      return;
    }

    this.ctx!.save();
    this.ctx!.clearRect(x, y, w, h);
    this.ctx!.translate(x, y);
    this.ctx!.rect(0, 0, w, h);
    this.ctx!.clip();

    // Fixed dimensions
    const positionMarkersHeight = 30; // Fixed height for position markers
    const sliderHeight = this.config.sliderHeight;
    const minHeightForBasicHeader = positionMarkersHeight + sliderHeight + 5; // Basic header with markers and slider
    const minHeightForConservation = minHeightForBasicHeader + 40; // Additional height needed to show conservation

    // Calculate total available height
    const canvasWidth = w;
    const canvasHeight = h;

    // Calculate positions from bottom up
    const sliderTop = canvasHeight - sliderHeight;
    const positionMarkersTop = sliderTop - positionMarkersHeight;

    // Determine if we have enough height for conservation plot
    const hasSpaceForConservation = canvasHeight >= minHeightForConservation;
    const conservationHeight = hasSpaceForConservation ?
      this.config.conservationHeight || 40 : 0;

    // Only render conservation if we have enough space
    if (hasSpaceForConservation && this.config.conservationData && this.config.conservationData.length > 0) {
      // Calculate conservation top position - ensure it's directly above position markers
      const conservationTop = positionMarkersTop - conservationHeight;
      this.drawConservationPlot(0, conservationTop, canvasWidth, conservationHeight);
    }

    // Draw slider at the bottom
    this.ctx!.fillStyle = this.config.sliderColor;
    this.ctx!.fillRect(0, sliderTop, canvasWidth, sliderHeight);

    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = Math.max(1, this.config.windowStartPosition);
    const totalSliderRange = this.config.totalPositions - visiblePositionsN;

    const sliderStartPX = totalSliderRange <= 0 ? 0 :
      windowStart / (totalSliderRange) * (canvasWidth - this.sliderWidth);
    const sliderLengthPX = totalSliderRange <= 0 ? canvasWidth : this.sliderWidth;

    // Draw slider handle
    this.ctx!.fillStyle = 'rgba(150, 150, 150, 0.5)';
    this.ctx!.fillRect(sliderStartPX, sliderTop, sliderLengthPX, sliderHeight);

    // Draw position indicator on slider
    if (this.config.currentPosition >= 1 && this.config.currentPosition <= this.config.totalPositions) {
      const currentPositionRatio = (this.config.currentPosition - 1) / (this.config.totalPositions - 1);
      const notchX = Math.round(currentPositionRatio * canvasWidth);

      this.ctx!.fillStyle = '#3CB173';
      this.ctx!.fillRect(
        notchX - 1,
        sliderTop - 2,
        3,
        sliderHeight + 4
      );
    }

    // Draw position markers
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i;
      if (position > this.config.totalPositions) break;

      const x = i * this.config.positionWidth;
      const cellWidth = this.config.positionWidth;
      const cellCenterX = x + cellWidth / 2;

      // Draw cell background
      if (this.config.cellBackground) {
        this.ctx!.fillStyle = i % 2 === 0 ? 'rgba(248, 248, 248, 0.3)' : 'rgba(242, 242, 242, 0.2)';
        this.ctx!.fillRect(x, positionMarkersTop, cellWidth, positionMarkersHeight);

        // Cell borders
        this.ctx!.strokeStyle = 'rgba(220, 220, 220, 0.7)';
        this.ctx!.beginPath();
        this.ctx!.moveTo(x, positionMarkersTop);
        this.ctx!.lineTo(x, sliderTop);
        this.ctx!.stroke();
      }

      // Position dot
      this.ctx!.fillStyle = '#999999';
      this.ctx!.beginPath();
      this.ctx!.arc(cellCenterX, positionMarkersTop + 10, 1, 0, Math.PI * 2);
      this.ctx!.fill();

      // Position numbers
      if (position === this.config.currentPosition || ((position === 1 || position % 10 === 0) &&
        Math.abs(position - this.config.currentPosition) > 1)) {
        this.ctx!.fillStyle = '#333333';
        this.ctx!.font = '12px monospace';
        this.ctx!.textAlign = 'center';
        this.ctx!.textBaseline = 'middle';
        this.ctx!.fillText(position.toString(), cellCenterX, positionMarkersTop + 20);
      }

      // Highlight current position
      if (position === this.config.currentPosition) {
        this.ctx!.fillStyle = 'rgba(60, 177, 115, 0.2)';
        this.ctx!.fillRect(
          x,
          positionMarkersTop,
          cellWidth,
          positionMarkersHeight
        );
      }
    }

    this.ctx!.restore();

    // Update event element
    preventable.preventDefault();
    this.eventElement.style.display = 'block';
    this.eventElement.style.left = `${this.config.x}px`;
    this.eventElement.style.top = `${this.config.y}px`;
    this.eventElement.style.width = `${this.config.width}px`;
    this.eventElement.style.height = `${this.config.height}px`;
  }

  private getConservationColor(value: number): string {
    // Default color scheme: gradient from blue (low) to red (high)
    if (!this.config.conservationColorScheme || this.config.conservationColorScheme === 'default') {
      // Blue (low) -> green (medium) -> red (high)
      if (value < 0.33) {
        return `rgba(0, 0, 255, ${Math.max(0.3, value * 3)})`;
      } else if (value < 0.67) {
        return `rgba(0, ${Math.floor(255 * ((value - 0.33) * 3))}, ${Math.floor(255 * (1 - (value - 0.33) * 3))}, 0.8)`;
      } else {
        return `rgba(${Math.floor(255 * ((value - 0.67) * 3))}, ${Math.floor(255 * (1 - (value - 0.67) * 3))}, 0, 0.8)`;
      }
    }
    else if (this.config.conservationColorScheme === 'rainbow') {
      // Rainbow scheme: blue -> cyan -> green -> yellow -> red
      const hue = (1 - value) * 240; // 240 (blue) to 0 (red)
      return `hsla(${hue}, 100%, 50%, 0.8)`;
    }
    else if (this.config.conservationColorScheme === 'heatmap') {
      // Heatmap: black (low) -> red -> yellow -> white (high)
      if (value < 0.33) {
        const r = Math.floor(value * 3 * 255);
        return `rgba(${r}, 0, 0, 0.8)`;
      } else if (value < 0.67) {
        const g = Math.floor((value - 0.33) * 3 * 255);
        return `rgba(255, ${g}, 0, 0.8)`;
      } else {
        const b = Math.floor((value - 0.67) * 3 * 255);
        return `rgba(255, 255, ${b}, 0.8)`;
      }
    }

    // Fallback to grayscale
    return `rgba(0, 0, 0, ${value})`;
  }


  private drawConservationPlot(startX: number, startY: number, width: number, height: number): void {
    if (!this.ctx || !this.config.conservationData) return;

    const ctx = this.ctx;
    const visiblePositionsN = Math.floor(this.config.width / this.config.positionWidth);
    const windowStart = Math.max(1, this.config.windowStartPosition);

    // Draw background
    ctx.fillStyle = 'rgba(245, 245, 245, 0.3)';
    ctx.fillRect(startX, startY, width, height);

    // Draw title
    ctx.fillStyle = '#333333';
    ctx.font = '10px Arial';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    // ctx.fillText('Conservation', startX + 5, startY + 2);

    // Draw Y-axis and tick marks
    ctx.strokeStyle = '#999999';
    ctx.beginPath();
    ctx.moveTo(startX, startY + height - 1);
    ctx.lineTo(startX, startY + 12);
    ctx.stroke();

    ctx.fillStyle = '#999999';
    ctx.textAlign = 'right';
    ctx.font = '8px Arial';
    ctx.fillText('1.0', startX - 2, startY + 12);
    ctx.fillText('0.5', startX - 2, startY + height / 2);
    ctx.fillText('0.0', startX - 2, startY + height - 2);

    // Calculate bar width (slightly narrower than position width)
    const barWidth = Math.max(1, this.config.positionWidth - 2);

    // Draw the conservation bars
    for (let i = 0; i < visiblePositionsN; i++) {
      const position = windowStart + i - 1; // 0-based for array indexing
      if (position < 0 || position >= this.config.conservationData.length) continue;

      const conservationValue = this.config.conservationData[position];
      if (conservationValue === undefined || conservationValue === null) continue;

      // Calculate bar position and dimensions
      const barX = startX + (i * this.config.positionWidth) + 1;
      const barHeight = Math.max(1, (height - 12) * conservationValue);
      const barY = startY + height - barHeight;

      // Determine color based on conservation value
      const color = this.getConservationColor(conservationValue);

      // Draw the bar
      ctx.fillStyle = color;
      ctx.fillRect(barX, barY, barWidth, barHeight);

      // Add border for visibility
      if (barWidth > 3 && barHeight > 3) {
        ctx.strokeStyle = 'rgba(0, 0, 0, 0.2)';
        ctx.strokeRect(barX, barY, barWidth, barHeight);
      }
    }
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
