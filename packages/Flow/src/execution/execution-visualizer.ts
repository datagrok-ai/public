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
  [NodeExecStatus.idle]: {boxcolor: '#78909c', bgcolor: '#ffffff'},
  [NodeExecStatus.running]: {boxcolor: '#1976d2', bgcolor: '#e3f2fd'},   // Blue
  [NodeExecStatus.completed]: {boxcolor: '#43a047', bgcolor: '#ffffff'}, // Green
  [NodeExecStatus.errored]: {boxcolor: '#e53935', bgcolor: '#ffebee'},   // Red
  [NodeExecStatus.stale]: {boxcolor: '#9E9E9E', bgcolor: '#f5f5f5'},
};

const DOT_RADIUS = 5;

export class ExecutionVisualizer {
  private canvasController: CanvasController;
  private graphManager: GraphManager;
  private originalColors: Map<number, OriginalColors> = new Map();
  private nodeOverlays: Map<number, (ctx: CanvasRenderingContext2D) => void> = new Map();
  private runningCount: number = 0;

  constructor(canvasController: CanvasController, graphManager: GraphManager) {
    this.canvasController = canvasController;
    this.graphManager = graphManager;
  }

  highlightNode(nodeId: number, status: NodeExecStatus): void {
    const node = this.graphManager.graph.getNodeById(nodeId);
    if (!node) return;

    // Track running node count for animation loop
    const prevOverlay = this.nodeOverlays.get(nodeId);
    const wasRunning = prevOverlay !== undefined &&
      this.originalColors.has(nodeId);

    // Save original colors on first highlight
    if (!this.originalColors.has(nodeId)) {
      this.originalColors.set(nodeId, {
        color: node.color || '#546e7a',
        bgcolor: node.bgcolor || '#ffffff',
        boxcolor: node.boxcolor || '#78909c',
      });
    }

    const colors = STATUS_COLORS[status];
    node.boxcolor = colors.boxcolor;
    node.bgcolor = colors.bgcolor;

    // Install onDrawForeground overlay for status dot
    this.installOverlay(node, status);

    // Manage animation loop for running nodes
    if (status === NodeExecStatus.running) {
      this.runningCount++;
      this.canvasController.startAnimationLoop();
    } else if (wasRunning) {
      // Node transitioned away from running
    }

    // Check if any nodes are still running; if not, stop animation
    if (status !== NodeExecStatus.running)
      this.checkStopAnimation();

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
    this.runningCount = 0;
    this.canvasController.stopAnimationLoop();
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
    this.runningCount = 0;
    this.canvasController.stopAnimationLoop();
    this.canvasController.graphCanvas.setDirty(true, true);
  }

  private checkStopAnimation(): void {
    // Count actual running nodes
    let hasRunning = false;
    for (const [nodeId] of this.nodeOverlays) {
      const node = this.graphManager.graph.getNodeById(nodeId);
      if (node && node.boxcolor === STATUS_COLORS[NodeExecStatus.running].boxcolor) {
        hasRunning = true;
        break;
      }
    }
    if (!hasRunning) {
      this.runningCount = 0;
      this.canvasController.stopAnimationLoop();
    }
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
    // Position: same area as the collapse icon (top-left of title bar)
    const boxSize = 10;
    const x = boxSize * 0.5 + 1;
    const y = -(((LiteGraph as any).NODE_TITLE_HEIGHT || 30) - boxSize * 0.5 - 1);

    ctx.save();

    // For running state, pulse opacity
    if (status === NodeExecStatus.running) {
      const phase = (Date.now() % 1200) / 1200;
      ctx.globalAlpha = 0.4 + 0.6 * Math.sin(phase * Math.PI * 2);
    }

    ctx.beginPath();
    ctx.arc(x, y, DOT_RADIUS, 0, Math.PI * 2);
    ctx.fillStyle = STATUS_COLORS[status].boxcolor;
    ctx.fill();

    // Add a white border for contrast
    ctx.strokeStyle = 'rgba(255,255,255,0.8)';
    ctx.lineWidth = 1.5;
    ctx.stroke();

    // Error indicator: white "!" inside dot
    if (status === NodeExecStatus.errored) {
      ctx.fillStyle = '#fff';
      ctx.font = 'bold 8px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('!', x, y);
    }

    // Completed: white checkmark
    if (status === NodeExecStatus.completed) {
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(x - 2.5, y);
      ctx.lineTo(x - 0.5, y + 2.5);
      ctx.lineTo(x + 3, y - 2);
      ctx.stroke();
    }

    ctx.restore();
  }
}

// Need LiteGraph for NODE_TITLE_HEIGHT
import {LiteGraph} from 'litegraph.js';
