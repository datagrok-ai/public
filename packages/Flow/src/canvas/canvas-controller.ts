/* eslint-disable camelcase */
import {LiteGraph, LGraphCanvas, LGraphNode} from 'litegraph.js';
import {GraphManager} from './graph-manager';
import {areTypesCompatible} from '../types/type-map';
import {drawTitleBoxWithStatus} from '../execution/execution-visualizer';

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
  private tooltipEl!: HTMLDivElement;
  private _animationTimer: ReturnType<typeof setInterval> | null = null;

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
      autoresize: false, skip_render: true,
    });

    // Set Roboto font for all canvas text
    (this.graphCanvas as any).title_text_font =
      '' + LiteGraph.NODE_TEXT_SIZE + 'px \'Roboto\', \'Roboto Local\', Arial, sans-serif';
    (this.graphCanvas as any).inner_text_font =
      'normal ' + (LiteGraph as any).NODE_SUBTEXT_SIZE + 'px \'Roboto\', \'Roboto Local\', Arial, sans-serif';

    // Configure appearance - Spotfire-inspired light theme
    this.graphCanvas.render_curved_connections = true;
    this.graphCanvas.render_connection_arrows = true;
    this.graphCanvas.render_connections_border = false;
    this.graphCanvas.highquality_render = true;
    this.graphCanvas.render_shadows = true;
    this.graphCanvas.clear_background = true;
    this.graphCanvas.links_render_mode = 2; // SPLINE_LINK
    this.graphCanvas.allow_searchbox = false;
    (this.graphCanvas as any).show_info = false;
    (this.graphCanvas as any).render_canvas_border = false;

    // Spotfire-inspired: clean background with subtle dot grid
    this.graphCanvas.background_image = this.createDotGridImage();
    (this.graphCanvas as any).clear_background_color = '#ebedf2';
    // Prevent background darkening on zoom (LiteGraph overlays the tile on the solid fill)
    (this.graphCanvas as any).zoom_modify_alpha = false;
    (this.graphCanvas as any).editor_alpha = 1.0;
    this.graphCanvas.default_link_color = '#8892a0';
    (this.graphCanvas as any).connections_width = 2.5;

    // Dark title text — visible on both white (expanded) and colored (collapsed) backgrounds
    (this.graphCanvas as any).node_title_color = '#333';

    // Suppress the default LiteGraph selected box outline (we draw our own in patchNodeRendering)
    (this.graphCanvas as any).node_box_outline_color = 'transparent';

    // Patch node rendering for clean visuals + colored title bars
    this.patchNodeRendering();

    // Unified status/collapse circle: install onDrawTitleBox on all nodes
    this.patchTitleBox();

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

    // Canvas tooltip element for hover hints
    this.tooltipEl = document.createElement('div');
    this.tooltipEl.className = 'funcflow-canvas-tooltip';
    container.appendChild(this.tooltipEl);
    this.setupCanvasTooltips();

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

  /** Apply Spotfire-inspired light theme color constants globally */
  private applyLightThemeGlobals(): void {
    (LGraphCanvas as any).DEFAULT_BACKGROUND_COLOR = '#ebedf2';
    LiteGraph.NODE_DEFAULT_COLOR = '#BDBDBD';
    LiteGraph.NODE_DEFAULT_BGCOLOR = '#ffffff';
    LiteGraph.NODE_DEFAULT_BOXCOLOR = '#90a4ae';
    LiteGraph.NODE_TEXT_COLOR = '#333';
    LiteGraph.NODE_TITLE_COLOR = '#BDBDBD';
    (LiteGraph as any).NODE_SELECTED_TITLE_COLOR = '#1565c0';
    // Make the default box outline transparent — our patchNodeRendering handles selection border
    (LiteGraph as any).NODE_BOX_OUTLINE_COLOR = 'transparent';
    (LiteGraph as any).WIDGET_BGCOLOR = '#f5f5f5';
    (LiteGraph as any).WIDGET_TEXT_COLOR = '#333';
    (LiteGraph as any).WIDGET_OUTLINE_COLOR = 'transparent';
    (LiteGraph as any).LINK_COLOR = '#8892a0';
    (LiteGraph as any).EVENT_LINK_COLOR = '#8892a0';
    (LiteGraph as any).CONNECTING_LINK_COLOR = '#1976d2';
  }

  /** Replace LiteGraph's drawNodeShape for Spotfire-inspired visuals.
   *
   * For expanded nodes: draws colored title bar + white title text on top of original render.
   * For collapsed nodes: LiteGraph already draws a colored bar, we just set node_title_color='#333'
   * as safe fallback but the collapsed path in LiteGraph uses fgcolor (node.color) natively.
   * Also adds: soft shadow, separator fix, blue selection border. */
  private patchNodeRendering(): void {
    const origDrawNodeShape = (LGraphCanvas.prototype as any).drawNodeShape;
    (LGraphCanvas.prototype as any).drawNodeShape = function(
      node: any, ctx: CanvasRenderingContext2D, size: number[],
      fgcolor: string, bgcolor: string, selected: boolean, mouse_over: boolean,
    ) {
      const title_height = (LiteGraph as any).NODE_TITLE_HEIGHT;
      const rr = (this as any).round_radius || 10;
      const isCollapsed = node.flags?.collapsed;

      // 1) Draw soft shadow
      ctx.save();
      ctx.shadowColor = 'rgba(0,0,0,0.10)';
      ctx.shadowBlur = 10;
      ctx.shadowOffsetX = 0;
      ctx.shadowOffsetY = 3;
      ctx.beginPath();
      if (isCollapsed) {
        const cw = (node as any)._collapsed_width || size[0];
        (ctx as any).roundRect(0, -title_height, cw, title_height, [rr]);
      } else {
        (ctx as any).roundRect(0, -title_height, size[0], size[1] + title_height, [rr]);
      }
      ctx.fillStyle = 'rgba(255,255,255,0.01)';
      ctx.fill();
      ctx.restore();

      // 2) Call original (selected=false to suppress default gray outline)
      origDrawNodeShape.call(this, node, ctx, size, fgcolor, bgcolor, false, mouse_over);

      // 3) Fix separator line between title and body
      if (!isCollapsed) {
        ctx.fillStyle = bgcolor || '#ffffff';
        ctx.fillRect(0, -1, size[0] + 1, 2);
      }

      // 4) Selection border: blue when selected, subtle gray when not
      ctx.save();
      if (selected) {
        ctx.strokeStyle = '#1976d2';
        ctx.lineWidth = 1.5;
        ctx.globalAlpha = 1.0;
      } else {
        ctx.strokeStyle = '#b0bec5';
        ctx.lineWidth = 1;
        ctx.globalAlpha = 0.4;
      }
      ctx.beginPath();
      if (isCollapsed) {
        const cw = (node as any)._collapsed_width || size[0];
        (ctx as any).roundRect(0, -title_height, cw, title_height, [rr]);
      } else {
        (ctx as any).roundRect(0, -title_height, size[0], size[1] + title_height, [rr]);
      }
      ctx.stroke();
      ctx.restore();
    };
  }

  /** Install unified status/collapse circle on all nodes via LGraphNode prototype */
  private patchTitleBox(): void {
    (LGraphNode.prototype as any).onDrawTitleBox = function(
      ctx: CanvasRenderingContext2D, title_height: number, size: number[], scale: number,
    ) {
      drawTitleBoxWithStatus(ctx, this, title_height, size, scale);
    };
  }

  /** Detects hover over node collapse icon and I/O slots; shows cursor + tooltip */
  private setupCanvasTooltips(): void {
    const titleH = (LiteGraph as any).NODE_TITLE_HEIGHT || 30;
    const boxSize = 10;
    let lastTip = '';

    this.canvas.addEventListener('mousemove', (e: MouseEvent) => {
      const gc = this.graphCanvas as any;
      const ds = gc.ds;
      // Convert mouse to graph coordinates
      const rect = this.canvas.getBoundingClientRect();
      const mx = (e.clientX - rect.left) / ds.scale - ds.offset[0];
      const my = (e.clientY - rect.top) / ds.scale - ds.offset[1];

      let tip = '';
      let cursor = '';

      const graph = this.graphManager.graph;
      const nodes = (graph as any)._nodes as LGraphNode[];
      if (nodes) {
        for (let i = nodes.length - 1; i >= 0; i--) {
          const node = nodes[i];
          const nx = node.pos[0];
          const ny = node.pos[1];

          // Check collapse icon area (top-left of title bar)
          const iconX = nx + boxSize * 0.5;
          const iconY = ny - titleH + boxSize * 0.5;
          if (mx >= iconX - 6 && mx <= iconX + 12 && my >= iconY - 6 && my <= iconY + 12) {
            tip = node.flags?.collapsed ? 'Expand node' : 'Collapse node';
            cursor = 'pointer';
            break;
          }

          // Account for collapsed nodes — they only show the title bar
          const isCollapsed = !!(node.flags as any)?.collapsed;
          const nodeW = isCollapsed ?
            ((node as any)._collapsed_width || node.size[0]) :
            node.size[0];

          // Check if over node title bar (for general info)
          if (mx >= nx && mx <= nx + nodeW && my >= ny - titleH && my < ny) {
            // Over title — show node description if it's a func node
            const desc = (node as any).dgFunc?.description;
            if (desc) {
              tip = desc;
              break;
            }
          }

          // Skip I/O slot tooltips when collapsed (slots are not visible)
          if (!isCollapsed) {
            // Check I/O slots
            if (node.inputs) {
              for (let s = 0; s < node.inputs.length; s++) {
                const slotY = ny + (LiteGraph as any).NODE_TITLE_HEIGHT * 0.5 + s * (LiteGraph as any).NODE_SLOT_HEIGHT;
                if (Math.abs(mx - nx) < 12 && Math.abs(my - slotY) < 8) {
                  const inp = node.inputs[s];
                  tip = `${inp.name} (${inp.type})`;
                  break;
                }
              }
            }
            if (!tip && node.outputs) {
              const ptCount = (node as any).properties?.['_passthroughCount'] ?? 0;
              for (let s = 0; s < node.outputs.length; s++) {
                const slotY = ny + (LiteGraph as any).NODE_TITLE_HEIGHT * 0.5 +
                  s * (LiteGraph as any).NODE_SLOT_HEIGHT;
                if (Math.abs(mx - (nx + node.size[0])) < 12 && Math.abs(my - slotY) < 8) {
                  const out = node.outputs[s];
                  if (s < ptCount) {
                    const inputName = node.inputs?.[s]?.name || out.name;
                    tip = `${inputName} \u2014 pass-through (${out.type}). Connect to enforce execution order`;
                  } else
                    tip = `${out.name} (${out.type})`;
                  break;
                }
              }
            }
          }

          if (tip) break;
        }
      }

      // Update cursor
      if (cursor)
        this.canvas.style.cursor = cursor || '';

      // Update tooltip
      if (tip !== lastTip) {
        lastTip = tip;
        if (tip) {
          this.tooltipEl.textContent = tip;
          this.tooltipEl.style.display = 'block';
        } else
          this.tooltipEl.style.display = 'none';
      }
      if (tip) {
        this.tooltipEl.style.left = (e.clientX - rect.left + 12) + 'px';
        this.tooltipEl.style.top = (e.clientY - rect.top + 12) + 'px';
      }
    });

    this.canvas.addEventListener('mouseleave', () => {
      this.tooltipEl.style.display = 'none';
      this.canvas.style.cursor = '';
      lastTip = '';
    });
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

  /** Start a periodic dirty flag to force canvas redraws for animations */
  startAnimationLoop(): void {
    if (this._animationTimer) return;
    this._animationTimer = setInterval(() => {
      this.graphCanvas.setDirty(true, true);
    }, 60); // ~16fps for smooth pulsing
  }

  /** Stop the periodic animation redraw */
  stopAnimationLoop(): void {
    if (this._animationTimer) {
      clearInterval(this._animationTimer);
      this._animationTimer = null;
    }
  }

  /** Generate a small tiled dot-grid image as a data URL for the canvas background */
  private createDotGridImage(): string {
    const size = 20;
    const sq = 1; // half-size of the square
    const cx = size / 2;
    const cy = size / 2;
    const c = document.createElement('canvas');
    c.width = size;
    c.height = size;
    const ctx = c.getContext('2d')!;
    ctx.fillStyle = '#ebedf2';
    ctx.fillRect(0, 0, size, size);
    ctx.fillStyle = '#8b8b8b';
    ctx.fillRect(cx - sq, cy - sq, sq * 2, sq * 2);
    return c.toDataURL();
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
