
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {TreeNode, Tree} from './aizynth-api';
import {HORIZONTAL_SPACING, MOL_HEIGHT, MOL_WIDTH, PADDING, VERTICAL_SPACING} from './const';

export function createPathsTreeTabs(paths: Tree[], isHorizontal: boolean = true): DG.TabControl {
  const tabControl = ui.tabControl();

  for (let i = 0; i < paths.length; i++) {
    tabControl.addPane(paths[i].scores['state score'].toFixed(3),
      () => {
        const pathObject = buildNestedStructure(paths[i]);
        const container = ui.div([], {classes: 'retrosynthesis-reaction-container'});
        createReactionTree(pathObject, isHorizontal).then((c) => {
          container.appendChild(c);
          //this.rotate !== 'None' && this.rotateCanvas90Degrees(c, this.rotate === 'Clockwise');
        });
        return container;
      });
  }
  return tabControl;
}


export function buildNestedStructure(node: TreeNode | Tree): { [key: string]: any } {
  const result: { [key: string]: any } = {};

  // If the current node is of type 'mol', add it to the result
  if (node.type === 'mol') {
    const smiles = node.smiles;
    const children = node.children || [];

    // Recursively process children
    const nestedChildren: { [key: string]: any } = {};
    children.forEach((child) => {
      const childResult = buildNestedStructure(child);
      Object.assign(nestedChildren, childResult);
    });

    // Add the current node's smiles as a key and its nested children as the value
    result[smiles] = nestedChildren;
  } else if (node.type === 'reaction') {
    // If the current node is of type 'reaction', process its children
    const children = node.children || [];
    children.forEach((child) => {
      const childResult = buildNestedStructure(child);
      Object.assign(result, childResult);
    });
  }
  return result;
}


export async function createReactionTree(reactionData: any, isHorizontal: boolean = true): Promise<HTMLCanvasElement> {
  const moleculesPerLevel: any = {};

  // Traverse the reaction data to determine the number of molecules per level
  function traverseLevels(data: any, level = 0) {
    Object.keys(data).forEach((product) => {
      moleculesPerLevel[level] ??= [];
      moleculesPerLevel[level].push(product);
      traverseLevels(data[product], level + 1);
    });
  }
  traverseLevels(reactionData);

  const maxLevel = Object.keys(moleculesPerLevel).length - 1;
  const levelHeights = new Array(maxLevel + 1).fill(0);
  const levelWidths = new Array(maxLevel + 1).fill(0);

  // Calculate the height and width for each level
  Object.entries(moleculesPerLevel).forEach(([level, molecules]) => {
    const numMolecules = (molecules as string[]).length;
    if (isHorizontal) {
      levelHeights[Number(level)] =
        numMolecules * (MOL_HEIGHT + PADDING) + (numMolecules - 1) * VERTICAL_SPACING + PADDING * 2;
      levelWidths[Number(level)] = MOL_WIDTH + PADDING * 2;
    } else {
      levelWidths[Number(level)] =
        numMolecules * (MOL_WIDTH + PADDING) + (numMolecules - 1) * HORIZONTAL_SPACING + PADDING * 2;
      levelHeights[Number(level)] = MOL_HEIGHT + PADDING * 2;
    }
  });

  // Calculate the total canvas size
  const neededHeight = isHorizontal ?
    levelHeights.reduce((acc, height) => Math.max(acc, height), 0) + PADDING * 2 :
    levelHeights.reduce((acc, height) => acc + height + VERTICAL_SPACING, 0);
  const neededWidth = isHorizontal ?
    levelWidths.reduce((acc, width) => acc + width + HORIZONTAL_SPACING, 0) :
    levelWidths.reduce((acc, width) => Math.max(acc, width), 0) + PADDING * 2;

  const canvas = document.createElement('canvas');
  canvas.width = neededWidth;
  canvas.height = neededHeight;
  const ctx = canvas.getContext('2d')!;
  const hitPositions: any[] = [];
  const reagentStartHeights = new Array(maxLevel + 1).fill(-1);
  const reagentStartWidths = new Array(maxLevel + 1).fill(-1);

  await drawTree(ctx, levelHeights, levelWidths, reagentStartHeights, reagentStartWidths,
    neededHeight, neededWidth, hitPositions, reactionData, 0, 0, 0, null, isHorizontal);
  return canvas;
}

async function drawTree(
  ctx: CanvasRenderingContext2D,
  levelHeights: number[],
  levelWidths: number[],
  reagentStartHeights: number[],
  reagentStartWidths: number[],
  neededHeight: number,
  neededWidth: number,
  hitPositions: any[],
  data: any,
  startX: number,
  startY: number,
  level: number,
  parent: string | null,
  isHorizontal: boolean,
) {
  const color = parent ? stringToColor(parent) : '#2083D5';
  const molecules = Object.keys(data);

  const levelStartX = isHorizontal ?
    level * (MOL_WIDTH + HORIZONTAL_SPACING + PADDING) + PADDING :
    (neededWidth - levelWidths[level]) / 2;
  const levelStartY = isHorizontal ?
    (neededHeight - levelHeights[level]) / 2 :
    level * (MOL_HEIGHT + VERTICAL_SPACING + PADDING) + PADDING;

  if (reagentStartHeights[level] === -1) reagentStartHeights[level] = levelStartY;
  if (reagentStartWidths[level] === -1) reagentStartWidths[level] = levelStartX;

  for (const mol of molecules) {
    if (level > 0) {
      if (isHorizontal)
        reactionLinesToMolecule(ctx, startX, startY, levelStartX, reagentStartHeights[level], isHorizontal);
      else
        reactionLinesToMolecule(ctx, startX, startY, reagentStartWidths[level], levelStartY, isHorizontal);
    }

    const {width, height} = await drawMolecule(ctx, mol, reagentStartWidths[level], reagentStartHeights[level], color);
    hitPositions.push({x: reagentStartWidths[level], y: reagentStartHeights[level], width, height, mol});

    await drawTree(
      ctx,
      levelHeights,
      levelWidths,
      reagentStartHeights,
      reagentStartWidths,
      neededHeight,
      neededWidth,
      hitPositions,
      data[mol],
      isHorizontal ? levelStartX : reagentStartWidths[level],
      isHorizontal ? reagentStartHeights[level] : levelStartY,
      level + 1,
      mol,
      isHorizontal,
    );

    if (isHorizontal)
      reagentStartHeights[level] += height + VERTICAL_SPACING;
    else
      reagentStartWidths[level] += width + HORIZONTAL_SPACING;
  }
}

function reactionLinesToMolecule(
  ctx: CanvasRenderingContext2D,
  startX: number,
  startY: number,
  endX: number,
  endY: number,
  isHorizontal: boolean,
) {
  if (isHorizontal) {
    // Horizontal orientation: draw lines from right to left
    const startLineX = startX + MOL_WIDTH + PADDING / 2;
    const startLineY = startY + MOL_HEIGHT / 2;
    const endLineX = endX - PADDING / 2;
    const endLineY = endY + MOL_HEIGHT / 2;
    const midX = (startLineX + endLineX) / 2;

    drawLine(ctx, startLineX, startLineY, midX, startLineY, true);
    drawLine(ctx, midX, startLineY, midX, endLineY);
    drawLine(ctx, midX, endLineY, endLineX, endLineY);
  } else {
    // Vertical orientation: draw lines from top to bottom
    const startLineX = startX + MOL_WIDTH / 2;
    const startLineY = startY + MOL_HEIGHT + PADDING / 2;
    const endLineX = endX + MOL_WIDTH / 2;
    const endLineY = endY - PADDING / 2;
    const midY = (startLineY + endLineY) / 2;

    drawLine(ctx, startLineX, startLineY, startLineX, midY, true);
    drawLine(ctx, startLineX, midY, endLineX, midY);
    drawLine(ctx, endLineX, midY, endLineX, endLineY);
  }
}

// Helper functions (unchanged)
function stringToColor(str: string): string {
  let hash = 0;
  for (let i = 0; i < str.length; i++) hash = str.charCodeAt(i) + ((hash << 5) - hash);
  let color = '#';
  for (let i = 0; i < 3; i++) {
    const value = (hash >> (i * 8)) & 0xff;
    color += ('00' + value.toString(16)).substr(-2);
  }
  return color;
}

function drawLine(ctx: CanvasRenderingContext2D, startX: number, startY: number, endX: number,
  endY: number, withDot = false) {
  ctx.save();
  ctx.beginPath();
  ctx.strokeStyle = 'black';
  ctx.moveTo(startX, startY);
  ctx.lineTo(endX, endY);
  if (withDot) ctx.arc(endX, endY, 2, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();
}

async function drawMolecule(ctx: CanvasRenderingContext2D, smiles: string, x: number, y: number, color: string) {
  const tempCanvas = document.createElement('canvas');
  tempCanvas.width = MOL_WIDTH;
  tempCanvas.height = MOL_HEIGHT;
  await grok.chem.canvasMol(0, 0, MOL_WIDTH, MOL_HEIGHT, tempCanvas, smiles, null, { //@ts-ignore
    autoCrop: true,
    autoCropMargin: 0,
    suppressChiralText: true,
  });
  ctx.drawImage(tempCanvas, x, y, MOL_WIDTH, MOL_HEIGHT);
  ctx.save();
  ctx.strokeStyle = color;
  ctx.beginPath();
  ctx.roundRect(x - PADDING / 2, y - PADDING / 2, MOL_WIDTH + PADDING, MOL_HEIGHT + PADDING, 5);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();
  return {width: MOL_WIDTH + PADDING, height: MOL_HEIGHT + PADDING};
}

// function rotateCanvas90Degrees(canvas: HTMLCanvasElement, clockwise = true) {
//   const ctx = canvas.getContext('2d')!;

//   // Create a new canvas
//   const newCanvas = ui.canvas(canvas.height, canvas.width);
//   const newCtx = newCanvas.getContext('2d')!;

//   // Set dimensions of new canvas
//   newCanvas.width = canvas.height;
//   newCanvas.height = canvas.width;

//   // Translate and rotate the new canvas
//   if (clockwise) {
//     newCtx.translate(canvas.height, 0);
//     newCtx.rotate(Math.PI / 2);
//   } else {
//     newCtx.translate(0, canvas.width);
//     newCtx.rotate(-Math.PI / 2);
//   }

//   // Draw the original canvas onto the new one
//   newCtx.drawImage(canvas, 0, 0);

//   // Clear the original canvas
//   ctx.clearRect(0, 0, canvas.width, canvas.height);

//   // Set new dimensions for the original canvas
//   canvas.width = newCanvas.width;
//   canvas.height = newCanvas.height;
//   canvas.style.width = `${newCanvas.width}px`;
//   canvas.style.height = `${newCanvas.height}px`;

//   // Draw the rotated image back to the original canvas
//   ctx.drawImage(newCanvas, 0, 0);
// }
