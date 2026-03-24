/** Maps execution state to LiteGraph node visual properties.
 *
 * The status indicator is unified with the LiteGraph "title box" circle
 * (top-left of every node, also used as collapse/expand toggle).
 *
 * Idle nodes: white fill + gray outline (drawn by patchTitleBox in CanvasController).
 * Status nodes: colored fill drawn via onDrawTitleBox override, with optional
 * icon (checkmark for completed, "!" for errored, pulse for running). */
import {LGraphNode, LiteGraph} from 'litegraph.js';
import {CanvasController} from '../canvas/canvas-controller';
import {GraphManager} from '../canvas/graph-manager';
import {NodeExecStatus} from './execution-state';

interface OriginalColors {
  color: string;
  bgcolor: string;
  boxcolor: string;
}

export const STATUS_COLORS: Record<NodeExecStatus, {boxcolor: string; bgcolor: string}> = {
  [NodeExecStatus.idle]: {boxcolor: '#90a4ae', bgcolor: '#ffffff'},
  [NodeExecStatus.running]: {boxcolor: '#1976d2', bgcolor: '#e3f2fd'},
  [NodeExecStatus.completed]: {boxcolor: '#43a047', bgcolor: '#ffffff'},
  [NodeExecStatus.errored]: {boxcolor: '#e53935', bgcolor: '#ffebee'},
  [NodeExecStatus.stale]: {boxcolor: '#9E9E9E', bgcolor: '#f5f5f5'},
};

/** Stored on each node so the title-box patch can read it */
const NODE_STATUS_KEY = '__ffStatus';

export class ExecutionVisualizer {
  private canvasController: CanvasController;
  private graphManager: GraphManager;
  private originalColors: Map<number, OriginalColors> = new Map();
  private trackedNodes: Set<number> = new Set();
  private runningCount: number = 0;

  constructor(canvasController: CanvasController, graphManager: GraphManager) {
    this.canvasController = canvasController;
    this.graphManager = graphManager;
  }

  highlightNode(nodeId: number, status: NodeExecStatus): void {
    const node = this.graphManager.graph.getNodeById(nodeId);
    if (!node) return;

    const wasRunning = (node as any)[NODE_STATUS_KEY] === NodeExecStatus.running;

    // Save original colors on first highlight
    if (!this.originalColors.has(nodeId)) {
      this.originalColors.set(nodeId, {
        color: node.color || '#BDBDBD',
        bgcolor: node.bgcolor || '#ffffff',
        boxcolor: node.boxcolor || '#90a4ae',
      });
    }

    const colors = STATUS_COLORS[status];
    node.boxcolor = colors.boxcolor;
    node.bgcolor = colors.bgcolor;
    (node as any)[NODE_STATUS_KEY] = status;
    this.trackedNodes.add(nodeId);

    // Manage animation loop for running nodes
    if (status === NodeExecStatus.running) {
      this.runningCount++;
      this.canvasController.startAnimationLoop();
    } else if (wasRunning) {
      // transitioned away
    }

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
        delete (node as any)[NODE_STATUS_KEY];
      }
    }
    this.originalColors.clear();
    this.trackedNodes.clear();
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
        (node as any)[NODE_STATUS_KEY] = NodeExecStatus.stale;
      }
    }
    this.runningCount = 0;
    this.canvasController.stopAnimationLoop();
    this.canvasController.graphCanvas.setDirty(true, true);
  }

  private checkStopAnimation(): void {
    let hasRunning = false;
    for (const nodeId of this.trackedNodes) {
      const node = this.graphManager.graph.getNodeById(nodeId);
      if (node && (node as any)[NODE_STATUS_KEY] === NodeExecStatus.running) {
        hasRunning = true;
        break;
      }
    }
    if (!hasRunning) {
      this.runningCount = 0;
      this.canvasController.stopAnimationLoop();
    }
  }
}

/** Draw the unified status/collapse circle.  Called from CanvasController.patchTitleBox(). */
export function drawTitleBoxWithStatus(
  ctx: CanvasRenderingContext2D, node: any, title_height: number, _size: number[], _scale: number,
): void {
  const status: NodeExecStatus | undefined = node[NODE_STATUS_KEY];
  const boxRadius = 5;
  const cx = title_height * 0.5;
  const cy = title_height * -0.5;

  ctx.save();

  if (status === NodeExecStatus.running) {
    const phase = (Date.now() % 1200) / 1200;
    ctx.globalAlpha = 0.4 + 0.6 * Math.sin(phase * Math.PI * 2);
  }

  // Fill circle
  ctx.beginPath();
  ctx.arc(cx, cy, boxRadius, 0, Math.PI * 2);
  if (!status || status === NodeExecStatus.idle) {
    // Idle: white fill + gray outline
    ctx.fillStyle = '#fff';
    ctx.fill();
    ctx.strokeStyle = '#90a4ae';
    ctx.lineWidth = 1.5;
    ctx.stroke();
  } else {
    ctx.fillStyle = STATUS_COLORS[status].boxcolor;
    ctx.fill();
  }

  // Status icons inside the circle
  if (status === NodeExecStatus.errored) {
    ctx.fillStyle = '#fff';
    ctx.font = 'bold 8px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('!', cx, cy);
  } else if (status === NodeExecStatus.completed) {
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(cx - 2.5, cy);
    ctx.lineTo(cx - 0.5, cy + 2.5);
    ctx.lineTo(cx + 3, cy - 2);
    ctx.stroke();
  }

  ctx.restore();
}
