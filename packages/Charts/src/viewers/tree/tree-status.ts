import * as DG from 'datagrok-api/dg';
import type {TreeViewer} from './tree-viewer';
import {Box, Readings, drawnItems, treePath} from '../../utils/utils';

/** The side of the square a branch's hit area reports, centred on a point of the branch line. */
const HIT = 4;

/** Where along a branch line, from its parent (0) to its node (1), a click target is looked for —
 * the middle first, away from both nodes' labels. */
const ALONG = [0.5, 0.35, 0.65, 0.2, 0.8];

type Branch = {path: string, value: number, point: {x: number, y: number} | null};

/** The line into a node is the click target: the viewer resolves a click to the node the clicked line
 * leads into, while a click on the node's symbol only folds or unfolds it. The point is kept only
 * where zrender's own hit test finds that line, not a label or another line on top. */
function branchPoint(edge: any, zr: any): {x: number, y: number} | null {
  const m = edge?.getComputedTransform?.();
  for (const t of ALONG) {
    const local = edge?.pointAt?.(t);
    if (local == null)
      return null;
    const x = m ? m[0] * local[0] + m[2] * local[1] + m[4] : local[0];
    const y = m ? m[1] * local[0] + m[3] * local[1] + m[5] : local[1];
    const hover = zr?.handler?.findHover?.(x, y);
    if (hover === undefined ? edge.contain?.(x, y) !== false : hover.target === edge)
      return {x, y};
  }
  return null;
}

/** The nodes below the root the layout drew — a folded node's descendants have no layout and no
 * element — addressed from the root label down: `All | false | F | Asian`. */
function drawnBranches(chart: any): Branch[] {
  const data = chart?.getModel?.()?.getSeriesByIndex?.(0)?.getData?.();
  const zr = chart?.getZr?.();
  const branches: Branch[] = [];
  for (const i of drawnItems(data)) {
    const node = data.tree?.getNodeByDataIndex?.(i);
    const names = treePath(node);
    if (names.length < 2 || data.getItemLayout(i) == null)
      continue;
    const point = branchPoint(data.getItemGraphicEl(i).__edge, zr);
    branches.push({path: names.join(' | '), value: Number(node.getValue?.()), point});
  }
  return branches;
}

/** While the viewer shows its message nothing of the previous frame is reported. */
export function treeStatus(v: TreeViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.chart.getDom().querySelector('canvas') : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['rows shown'] = v.filter.trueCount;
  const branches = canvas ? drawnBranches(v.chart) : [];
  values['branches'] = branches.length;
  for (const branch of branches) {
    values[`rows of branch ${branch.path}`] = branch.value;
    if (branch.point !== null) {
      hitAreas[`branch ${branch.path}`] = {
        x: branch.point.x - HIT / 2, y: branch.point.y - HIT / 2, width: HIT, height: HIT,
      };
    }
  }

  if (canvas) {
    const box = canvas.getBoundingClientRect();
    hitAreas['view'] = {x: 0, y: 0, width: box.width, height: box.height};
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
