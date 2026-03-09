/* eslint-disable camelcase */
import {LiteGraph, LGraphCanvas, LGraphNode} from 'litegraph.js';
import {GraphManager} from './graph-manager';
import {areTypesCompatible} from '../types/type-map';

export interface CanvasCallbacks {
  onNodeSelected?: (node: LGraphNode) => void;
  onNodeDeselected?: (node: LGraphNode) => void;
  onGraphChanged?: () => void;
}

/** Controls the LGraphCanvas: rendering, events, resize */
export class CanvasController {
  graphCanvas: LGraphCanvas;
  canvas: HTMLCanvasElement;
  graphManager: GraphManager;

  constructor(container: HTMLElement, graphManager: GraphManager, callbacks?: CanvasCallbacks) {
    this.graphManager = graphManager;

    // Override connection type validation BEFORE creating the canvas
    this.setupTypeValidation();

    // Create the canvas element
    this.canvas = document.createElement('canvas');
    this.canvas.style.width = '100%';
    this.canvas.style.height = '100%';
    this.canvas.style.display = 'block';
    container.appendChild(this.canvas);

    this.resizeToContainer(container);

    // Apply light theme globals BEFORE creating the canvas
    this.applyLightThemeGlobals();

    // Create LGraphCanvas
    this.graphCanvas = new LGraphCanvas(this.canvas, graphManager.graph, {
      autoresize: false,
    });

    // Configure appearance - light theme
    this.graphCanvas.render_curved_connections = true;
    this.graphCanvas.render_connection_arrows = false;
    this.graphCanvas.render_connections_border = false;
    this.graphCanvas.highquality_render = true;
    this.graphCanvas.render_shadows = false;
    this.graphCanvas.clear_background = true;
    this.graphCanvas.links_render_mode = 2; // SPLINE_LINK
    this.graphCanvas.allow_searchbox = true;
    (this.graphCanvas as any).show_info = false;
    (this.graphCanvas as any).render_canvas_border = false;

    // Light theme: the key properties for background rendering
    this.graphCanvas.background_image = '';
    (this.graphCanvas as any).clear_background_color = '#e8e8e8';
    this.graphCanvas.default_link_color = '#555';
    (this.graphCanvas as any).connections_width = 2;

    // Dark title text — readable on both white and colored title bars
    (this.graphCanvas as any).node_title_color = '#333';

    // Setup callbacks
    if (callbacks?.onNodeSelected)
      this.graphCanvas.onNodeSelected = callbacks.onNodeSelected;
    if (callbacks?.onNodeDeselected)
      this.graphCanvas.onNodeDeselected = callbacks.onNodeDeselected;
    if (callbacks?.onGraphChanged) {
      const graph = graphManager.graph;
      const onChange = callbacks.onGraphChanged;
      graph.onNodeAdded = () => onChange();
      graph.onNodeRemoved = () => onChange();
      graph.onLinkAdded = () => onChange();
      graph.onLinkRemoved = () => onChange();
    }

    // Observe container resize
    const resizeObserver = new ResizeObserver(() => {
      this.resizeToContainer(container);
      this.graphCanvas.setDirty(true, true);
    });
    resizeObserver.observe(container);
  }

  /** Override LiteGraph's type validation to support our custom types */
  private setupTypeValidation(): void {
    const lg = LiteGraph as any;
    lg.isValidConnection = function(type_a: string, type_b: string): boolean {
      // Let LiteGraph handle its native wildcard types (0, "", "*", null)
      if (!type_a || !type_b || type_a === '*' || type_b === '*')
        return true;
      if (type_a === type_b) return true;

      // Our custom "dynamic" type connects to anything
      if (type_a === 'dynamic' || type_b === 'dynamic') return true;
      if (type_a === 'object' || type_b === 'object') return true;

      // Use our type compatibility rules
      return areTypesCompatible(String(type_a), String(type_b));
    };
  }

  /** Apply light theme color constants globally */
  private applyLightThemeGlobals(): void {
    (LGraphCanvas as any).DEFAULT_BACKGROUND_COLOR = '#e8e8e8';
    LiteGraph.NODE_DEFAULT_COLOR = '#BDBDBD';
    LiteGraph.NODE_DEFAULT_BGCOLOR = '#ffffff';
    LiteGraph.NODE_DEFAULT_BOXCOLOR = 'transparent';
    LiteGraph.NODE_TEXT_COLOR = '#333';
    LiteGraph.NODE_TITLE_COLOR = '#BDBDBD';
    (LiteGraph as any).NODE_SELECTED_TITLE_COLOR = '#000';
    (LiteGraph as any).NODE_BOX_OUTLINE_COLOR = '#666';
    (LiteGraph as any).WIDGET_BGCOLOR = '#f5f5f5';
    (LiteGraph as any).WIDGET_TEXT_COLOR = '#333';
    (LiteGraph as any).WIDGET_OUTLINE_COLOR = 'transparent';
    (LiteGraph as any).LINK_COLOR = '#666';
    (LiteGraph as any).EVENT_LINK_COLOR = '#666';
    (LiteGraph as any).CONNECTING_LINK_COLOR = '#444';

    // this.patchNodeRendering();
  }

  /** Replace LiteGraph's drawNodeShape to fix light-theme rendering artifacts.
   *
   * The original method draws a semi-transparent separator line at the
   * title/body junction and leaves the node body without a stroke, which
   * produces faint white edge artifacts on a light background. This
   * replacement keeps all the original rendering but:
   *  - Skips the separator fillRect
   *  - Adds a subtle border stroke around the whole node shape */
  private patchNodeRendering(): void {
    const origDrawNodeShape = (LGraphCanvas.prototype as any).drawNodeShape;
    (LGraphCanvas.prototype as any).drawNodeShape = function(
      node: any, ctx: CanvasRenderingContext2D, size: number[],
      fgcolor: string, bgcolor: string, selected: boolean, mouse_over: boolean,
    ) {
      // Call the original rendering
      origDrawNodeShape.call(this, node, ctx, size, fgcolor, bgcolor, selected, mouse_over);

      // Overdraw the separator with the node bgcolor to suppress the artifact
      if (!node.flags.collapsed) {
        ctx.fillStyle = bgcolor || '#ffffff';
        ctx.fillRect(0, -1, size[0] + 1, 2);
      }

      // Draw a clean border around the node to mask any edge artifacts
      const title_height = (LiteGraph as any).NODE_TITLE_HEIGHT;
      const shape = node._shape || node.constructor.shape || (LiteGraph as any).ROUND_SHAPE;
      const rr = (this as any).round_radius || 10;

      ctx.save();
      ctx.globalAlpha = 0.25;
      ctx.strokeStyle = '#888';
      ctx.lineWidth = 1;
      ctx.beginPath();
      if (shape === (LiteGraph as any).BOX_SHAPE)
        ctx.rect(0, -title_height, size[0] + 1, size[1] + title_height);
      else
        (ctx as any).roundRect(0, -title_height, size[0] + 1, size[1] + title_height, [rr]);

      ctx.stroke();
      ctx.restore();
    };
  }

  private resizeToContainer(container: HTMLElement): void {
    const rect = container.getBoundingClientRect();
    if (rect.width > 0 && rect.height > 0) {
      this.canvas.width = rect.width;
      this.canvas.height = rect.height;
    }
  }

  zoomToFit(): void {
    this.graphCanvas.ds.reset();
    this.graphCanvas.setDirty(true, true);
  }

  zoomIn(): void {
    this.graphCanvas.ds.changeScale(1.2);
    this.graphCanvas.setDirty(true, true);
  }

  zoomOut(): void {
    this.graphCanvas.ds.changeScale(0.8);
    this.graphCanvas.setDirty(true, true);
  }

  startRendering(): void {
    this.graphCanvas.startRendering();
  }

  stopRendering(): void {
    this.graphCanvas.stopRendering();
  }

  addNodeAtCenter(nodeType: string): LGraphNode | null {
    try {
      const node = LiteGraph.createNode(nodeType);
      if (!node) return null;
      const area = this.graphCanvas.ds.visible_area;
      node.pos = [
        area[0] + area[2] / 2 - node.size[0] / 2,
        area[1] + area[3] / 2 - node.size[1] / 2,
      ];
      this.graphManager.graph.add(node);
      return node;
    } catch (e) {
      console.warn('Failed to add node:', e);
      return null;
    }
  }
}
