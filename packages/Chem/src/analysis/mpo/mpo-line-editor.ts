import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import Konva from 'konva';
import { PropertyDesirability } from './mpo';
import {Subject} from "rxjs"; // Import type from mpo.ts

// Constants for the editor layout
const EDITOR_PADDING = { top: 10, right: 10, bottom: 20, left: 30 };
const POINT_RADIUS = 3;


export class MpoDesirabilityLineEditor {
  root = ui.div();
  onChanged = new Subject();

  constructor(prop: PropertyDesirability, width: number, height: number) {
    this.root.style.width = `${width}px`;
    this.root.style.height = `${height}px`;
    this.root.style.position = 'relative'; // Needed for absolute positioning of Konva stage
    const that = this;

    // Delay Konva initialization slightly to ensure this.container is in DOM
    setTimeout(() => {
      if (!this.root.parentElement) {
        console.warn('Konva this.container not attached to DOM yet.');
        // Optionally, retry or handle error
        return;
      }

      const stage = new Konva.Stage({
        container: this.root, // Use the this.container div
        width: width,
        height: height,
      });

      const layer = new Konva.Layer();
      stage.add(layer);

      const minX = prop.min ?? Math.min(...prop.line.map(p => p[0]));
      const maxX = prop.max ?? Math.max(...prop.line.map(p => p[0]));

      // --- Draw Axes ---
      const xAxis = new Konva.Line({
        points: [EDITOR_PADDING.left, height - EDITOR_PADDING.bottom, width - EDITOR_PADDING.right, height - EDITOR_PADDING.bottom],
        stroke: 'grey', // Lighter color for axes
        strokeWidth: 1,
      });
      const yAxis = new Konva.Line({
        points: [EDITOR_PADDING.left, EDITOR_PADDING.top, EDITOR_PADDING.left, height - EDITOR_PADDING.bottom],
        stroke: 'grey', // Lighter color for axes
        strokeWidth: 1,
      });
      layer.add(xAxis, yAxis);

      // --- Draw Axis Labels/Ticks (Simplified) ---
      const minXLabel = new Konva.Text({
        x: EDITOR_PADDING.left,
        y: height - EDITOR_PADDING.bottom + 3, // Below axis
        text: minX.toFixed(1), // Fewer decimals for smaller space
        fontSize: 9,
        fill: 'grey'
      });
      const maxXLabel = new Konva.Text({
        x: width - EDITOR_PADDING.right - 15, // Adjust position
        y: height - EDITOR_PADDING.bottom + 3, // Below axis
        text: maxX.toFixed(1), // Fewer decimals
        fontSize: 9,
        align: 'right',
        fill: 'grey'
      });
      const zeroYLabel = new Konva.Text({
        x: EDITOR_PADDING.left - 20, // Left of axis
        y: height - EDITOR_PADDING.bottom - 5, // Align with axis bottom
        text: '0.0',
        fontSize: 9,
        fill: 'grey'
      });
      const oneYLabel = new Konva.Text({
        x: EDITOR_PADDING.left - 20, // Left of axis
        y: EDITOR_PADDING.top - 5, // Align with axis top
        text: '1.0',
        fontSize: 9,
        fill: 'grey'
      });
      layer.add(minXLabel, maxXLabel, zeroYLabel, oneYLabel);


      // --- Draw Line and Points ---
      const konvaLine = new Konva.Line({
        points: [], // Will be populated by points
        stroke: '#2077b4',
        strokeWidth: 2,
        lineCap: 'round',
        lineJoin: 'round'
      });
      layer.add(konvaLine);

      const pointsGroup = new Konva.Group(); // Group for draggable points
      layer.add(pointsGroup);

      // Function to redraw everything based on prop.line
      function redraw(notify: boolean = true) {
        pointsGroup.destroyChildren(); // Clear old points
        const konvaPoints: number[] = [];

        prop.line.sort((a, b) => a[0] - b[0]); // Ensure sorted

        prop.line.forEach((p, index) => {
          const coords = toCanvasCoords(p[0], p[1], minX, maxX, width, height);
          konvaPoints.push(coords.x, coords.y);

          const pointCircle = new Konva.Circle({
            x: coords.x,
            y: coords.y,
            radius: POINT_RADIUS,
            fill: '#d72f30',
            stroke: 'black',
            strokeWidth: 1,
            draggable: true, // Make points draggable
            hitStrokeWidth: 5 // Easier to hit for dragging/clicking
          });
          // Store index directly on the node for easy access
          pointCircle.setAttr('_pointIndex', index);

          // --- Dragging Logic ---
          pointCircle.on('dragmove', (evt: Konva.KonvaEventObject<DragEvent>) => {
            const circle = evt.target as Konva.Circle;
            const pos = circle.position();
            const currentPointIndex = circle.getAttr('_pointIndex');

            // Constrain dragging horizontally between neighbors (or bounds)
            const prevX = currentPointIndex > 0 ? prop.line[currentPointIndex - 1][0] : minX;
            const nextX = currentPointIndex < prop.line.length - 1 ? prop.line[currentPointIndex + 1][0] : maxX;
            // Add a small buffer to avoid points overlapping exactly
            const buffer = (maxX - minX === 0) ? 0 : 0.001 * (maxX - minX); // Avoid NaN if minX === maxX
            const minCanvasX = toCanvasCoords(prevX + (currentPointIndex > 0 ? buffer : 0), 0, minX, maxX, width, height).x;
            const maxCanvasX = toCanvasCoords(nextX - (currentPointIndex < prop.line.length - 1 ? buffer : 0), 0, minX, maxX, width, height).x;

            pos.x = Math.max(minCanvasX, Math.min(maxCanvasX, pos.x));

            // Constrain dragging vertically
            const plotTop = EDITOR_PADDING.top;
            const plotBottom = height - EDITOR_PADDING.bottom;
            pos.y = Math.max(plotTop, Math.min(plotBottom, pos.y));

            circle.position(pos); // Update position after constraints

            // Update data array
            const dataCoords = toDataCoords(pos.x, pos.y, minX, maxX, width, height);
            prop.line[currentPointIndex][0] = dataCoords.x;
            prop.line[currentPointIndex][1] = dataCoords.y;

            // Update the connecting line during drag
            const currentKonvaPoints = prop.line.map((pData, idx) => {
              // Use dragged circle position directly for the point being dragged
              if (idx === currentPointIndex) {
                return [pos.x, pos.y];
              } else {
                const c = toCanvasCoords(pData[0], pData[1], minX, maxX, width, height);
                return [c.x, c.y];
              }
            }).flat();
            konvaLine.points(currentKonvaPoints);
            layer.batchDraw(); // More efficient redraw
          });

          pointCircle.on('dragend', (evt: Konva.KonvaEventObject<DragEvent>) => {
            // Ensure data is sorted after drag, although constraints should handle it
            prop.line.sort((a, b) => a[0] - b[0]);
            redraw(); // Full redraw on drag end to fix indices and line path
          });

          // --- Right-click to Remove ---
          pointCircle.on('contextmenu', (evt: Konva.KonvaEventObject<MouseEvent>) => {
            evt.evt.preventDefault(); // Prevent browser context menu
            // Keep at least 2 points for a line segment
            if (prop.line.length <= 2) {
              // Use DG tooltip or simple alert/warning
              grok.shell.warning("Cannot remove points, minimum of 2 required.");
              return;
            }

            const circle = evt.target as Konva.Circle;
            const indexToRemove = circle.getAttr('_pointIndex');
            prop.line.splice(indexToRemove, 1);
            redraw(); // Redraw after removal
          });

          // Enhance usability: change cursor on hover and show tooltip
          pointCircle.on('mouseenter', (evt: Konva.KonvaEventObject<MouseEvent>) => {
            stage.container().style.cursor = 'pointer';
            const circle = evt.target as Konva.Circle;
            const pos = circle.position();
            const dataCoords = toDataCoords(pos.x, pos.y, minX, maxX, width, height);
            const tooltipText = `X: ${dataCoords.x.toFixed(2)}, Y: ${dataCoords.y.toFixed(2)}<br><br>Drag to move, right-click to delete`;
            ui.tooltip.show(tooltipText, evt.evt.clientX, evt.evt.clientY);
          });
          pointCircle.on('mouseleave', (evt: Konva.KonvaEventObject<MouseEvent>) => {
            stage.container().style.cursor = 'default';
            ui.tooltip.hide();
          });

          pointsGroup.add(pointCircle);
        });

        konvaLine.points(konvaPoints);
        layer.batchDraw();

        if (notify)
          that.onChanged.next();
      }

      // --- Left-click to Add Point ---
      stage.on('click tap', (evt: Konva.KonvaEventObject<PointerEvent>) => {
        // Ignore clicks on existing points (circles) or non-left clicks
        if (evt.target instanceof Konva.Circle || evt.evt.button !== 0) return;

        const pos = stage.getPointerPosition();
        if (!pos) return;

        // Ensure click is within the plot area boundaries
        if (pos.x < EDITOR_PADDING.left || pos.x > width - EDITOR_PADDING.right ||
          pos.y < EDITOR_PADDING.top || pos.y > height - EDITOR_PADDING.bottom) {
          return;
        }

        const dataCoords = toDataCoords(pos.x, pos.y, minX, maxX, width, height);
        // Add the new point
        prop.line.push([dataCoords.x, dataCoords.y]);

        // No need to sort here, redraw() handles sorting
        redraw();
      });

      // Change cursor when over the stage plot area for adding points
      stage.on('mouseenter', (evt: Konva.KonvaEventObject<MouseEvent>) => {
        if (!(evt.target instanceof Konva.Circle)) {
          const pos = stage.getPointerPosition();
          if (pos && pos.x >= EDITOR_PADDING.left && pos.x <= width - EDITOR_PADDING.right &&
            pos.y >= EDITOR_PADDING.top && pos.y <= height - EDITOR_PADDING.bottom) {
            stage.container().style.cursor = 'crosshair';
          }
        }
      });

      stage.on('mouseleave', (evt: Konva.KonvaEventObject<MouseEvent>) => {
        if (!(evt.target instanceof Konva.Circle)) {
          stage.container().style.cursor = 'default';
        }
      });

      stage.on('mousemove', (evt: Konva.KonvaEventObject<MouseEvent>) => {
        if (!(evt.target instanceof Konva.Circle)) {
          const pos = stage.getPointerPosition();
          if (pos && pos.x >= EDITOR_PADDING.left && pos.x <= width - EDITOR_PADDING.right &&
            pos.y >= EDITOR_PADDING.top && pos.y <= height - EDITOR_PADDING.bottom) {
            stage.container().style.cursor = 'crosshair';
          }
          else
            stage.container().style.cursor = 'default';
        }
      });

      stage.on('mouseout', (_) => ui.tooltip.hide());

      // Initial draw
      redraw(false);

    }, 0);
  }
}


// Function to transform data coordinates to canvas coordinates
function toCanvasCoords(x: number, y: number, minX: number, maxX: number, width: number, height: number): { x: number, y: number } {
  const plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
  const plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;
  // Handle case where minX === maxX to avoid division by zero
  const scaleX = (maxX - minX === 0) ? 1 : plotWidth / (maxX - minX);
  const scaleY = plotHeight; // y data is 0-1

  const canvasX = EDITOR_PADDING.left + (x - minX) * scaleX;
  const canvasY = EDITOR_PADDING.top + plotHeight - (y * scaleY); // Flip y-axis

  return { x: canvasX, y: canvasY };
}

// Function to transform canvas coordinates to data coordinates
function toDataCoords(canvasX: number, canvasY: number, minX: number, maxX: number, width: number, height: number): { x: number, y: number } {
  const plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
  const plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;
  // Handle case where minX === maxX
  const scaleX = (maxX - minX === 0) ? 1 : plotWidth / (maxX - minX);
  const scaleY = plotHeight;

  let dataX = minX + (canvasX - EDITOR_PADDING.left) / scaleX;
  let dataY = (EDITOR_PADDING.top + plotHeight - canvasY) / scaleY;

  // Clamp values
  dataX = Math.max(minX, Math.min(maxX, dataX));
  dataY = Math.max(0, Math.min(1, dataY));

  return { x: dataX, y: dataY };
}