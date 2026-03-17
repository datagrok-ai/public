/** Maps execution state to LiteGraph node visual properties */
import {LGraphNode} from 'litegraph.js';
import {CanvasController} from '../canvas/canvas-controller';
import {GraphManager} from '../canvas/graph-manager';
import {NodeExecStatus} from './execution-state';

interface OriginalColors {
  color: string;
  bgcolor: string;
  boxcolor: string;
}

const STATUS_COLORS: Record<NodeExecStatus, {boxcolor: string; bgcolor: string}> = {
  [NodeExecStatus.idle]: {boxcolor: '#888', bgcolor: '#ffffff'},
  [NodeExecStatus.running]: {boxcolor: '#FFA000', bgcolor: '#FFF8E1'},
  [NodeExecStatus.completed]: {boxcolor: '#4CAF50', bgcolor: '#ffffff'},
  [NodeExecStatus.errored]: {boxcolor: '#F44336', bgcolor: '#FFEBEE'},
  [NodeExecStatus.stale]: {boxcolor: '#9E9E9E', bgcolor: '#f5f5f5'},
};

const DOT_RADIUS = 5;
const DOT_MARGIN = 8;

export class ExecutionVisualizer {
  private canvasController: CanvasController;
  private graphManager: GraphManager;
  private originalColors: Map<number, OriginalColors> = new Map();
  private nodeOverlays: Map<number, (ctx: CanvasRenderingContext2D) => void> = new Map();

  constructor(canvasController: CanvasController, graphManager: GraphManager) {
    this.canvasController = canvasController;
    this.graphManager = graphManager;
  }

  highlightNode(nodeId: number, status: NodeExecStatus): void {
    const node = this.graphManager.graph.getNodeById(nodeId);
    if (!node) return;

    // Save original colors on first highlight
    if (!this.originalColors.has(nodeId)) {
      this.originalColors.set(nodeId, {
        color: node.color || '#BDBDBD',
        bgcolor: node.bgcolor || '#ffffff',
        boxcolor: node.boxcolor || '#888',
      });
    }

    const colors = STATUS_COLORS[status];
    node.boxcolor = colors.boxcolor;
    node.bgcolor = colors.bgcolor;

    // Install onDrawForeground overlay for status dot
    this.installOverlay(node, status);

    this.canvasController.graphCanvas.setDirty(true, true);
  }

  resetAllNodes(): void {
    for (const [nodeId, colors] of this.originalColors) {
      const node = this.graphManager.graph.getNodeById(nodeId);
      if (node) {
        node.color = colors.color;
        node.bgcolor = colors.bgcolor;
        node.boxcolor = colors.boxcolor;
        node.onDrawForeground = undefined as any;
      }
    }
    this.originalColors.clear();
    this.nodeOverlays.clear();
    this.canvasController.graphCanvas.setDirty(true, true);
  }

  markAllStale(): void {
    const staleColors = STATUS_COLORS[NodeExecStatus.stale];
    for (const [nodeId] of this.originalColors) {
      const node = this.graphManager.graph.getNodeById(nodeId);
      if (node) {
        node.boxcolor = staleColors.boxcolor;
        node.bgcolor = staleColors.bgcolor;
        this.installOverlay(node, NodeExecStatus.stale);
      }
    }
    this.canvasController.graphCanvas.setDirty(true, true);
  }

  private installOverlay(node: LGraphNode, status: NodeExecStatus): void {
    if (status === NodeExecStatus.idle) {
      node.onDrawForeground = undefined as any;
      this.nodeOverlays.delete(node.id);
      return;
    }

    const drawFn = (ctx: CanvasRenderingContext2D) => {
      this.drawStatusDot(ctx, node, status);
    };
    this.nodeOverlays.set(node.id, drawFn);
    node.onDrawForeground = drawFn;
  }

  private drawStatusDot(ctx: CanvasRenderingContext2D, node: LGraphNode, status: NodeExecStatus): void {
    // When collapsed, the node only renders as a title bar — use collapsed width
    const isCollapsed = !!(node.flags as any)?.collapsed;
    const nodeWidth = isCollapsed ?
      ((node as any)._collapsed_width || node.size[0]) :
      node.size[0];

    const x = nodeWidth - DOT_MARGIN;
    const y = -DOT_MARGIN;

    ctx.save();

    // For running state, pulse opacity
    if (status === NodeExecStatus.running) {
      const phase = (Date.now() % 1000) / 1000;
      ctx.globalAlpha = 0.5 + 0.5 * Math.sin(phase * Math.PI * 2);
    }

    ctx.beginPath();
    ctx.arc(x, y, DOT_RADIUS, 0, Math.PI * 2);
    ctx.fillStyle = STATUS_COLORS[status].boxcolor;
    ctx.fill();

    // Error indicator: white "!" inside dot
    if (status === NodeExecStatus.errored) {
      ctx.fillStyle = '#fff';
      ctx.font = 'bold 8px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('!', x, y);
    }

    ctx.restore();
  }
}
