
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TreeNode, Tree} from './aizynth-api';
import {HORIZONTAL_SPACING, MOL_HEIGHT, MOL_WIDTH, PADDING, VERTICAL_SPACING} from './const';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import '../css/aizynthfinder.css';

export function createPathsTreeTabs(paths: Tree[], isHorizontal: boolean = true): DG.TabControl {
  const tabControl = ui.tabControl();

  for (let i = 0; i < paths.length; i++) {
    tabControl.addPane(paths[i].scores['state score'].toFixed(3),
      () => {
        const pathObject = buildNestedStructure(paths[i]);
        const container = ui.div([], {classes: 'retrosynthesis-reaction-container'});
        createReactionTreeSVG(pathObject, isHorizontal).then((svgTree) => {
          container.appendChild(svgTree);
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


export function isFragment(molString: string) {
  if (DG.chem.isMolBlock(molString))
    return MolfileHandler.getInstance(molString).isFragment();
  else
    return !!molString.match(/\[.?:|\*.?\]/g);
}


export async function createReactionTreeSVG(reactionData: any, isHorizontal: boolean = true): Promise<SVGElement> {
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

  // Calculate the total SVG size
  const neededHeight = isHorizontal ?
    levelHeights.reduce((acc, height) => Math.max(acc, height), 0) + PADDING * 2 :
    levelHeights.reduce((acc, height) => acc + height + VERTICAL_SPACING, 0);
  const neededWidth = isHorizontal ?
    levelWidths.reduce((acc, width) => acc + width + HORIZONTAL_SPACING, 0) :
    levelWidths.reduce((acc, width) => Math.max(acc, width), 0) + PADDING * 2;

  // Create SVG element
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  setAttributes(svg, {'width': neededWidth.toString(), 'height': neededHeight.toString(),
    'viewBox': `0 0 ${neededWidth} ${neededHeight}`});

  // Add background
  const bgRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  setAttributes(bgRect, {'width': '100%', 'height': '100%', 'fill': 'white'});
  svg.appendChild(bgRect);

  const reagentStartHeights = new Array(maxLevel + 1).fill(-1);
  const reagentStartWidths = new Array(maxLevel + 1).fill(-1);

  await drawTreeSVG(svg, levelHeights, levelWidths, reagentStartHeights, reagentStartWidths,
    neededHeight, neededWidth, reactionData, 0, 0, 0, null, isHorizontal);

  return svg;
}

async function drawTreeSVG(
  svg: SVGElement,
  levelHeights: number[],
  levelWidths: number[],
  reagentStartHeights: number[],
  reagentStartWidths: number[],
  neededHeight: number,
  neededWidth: number,
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
        reactionLinesToMoleculeSVG(svg, startX, startY, levelStartX, reagentStartHeights[level], isHorizontal);
      else
        reactionLinesToMoleculeSVG(svg, startX, startY, reagentStartWidths[level], levelStartY, isHorizontal);
    }

    await drawMoleculeSVG(svg, mol, reagentStartWidths[level], reagentStartHeights[level], color);

    await drawTreeSVG(
      svg,
      levelHeights,
      levelWidths,
      reagentStartHeights,
      reagentStartWidths,
      neededHeight,
      neededWidth,
      data[mol],
      isHorizontal ? levelStartX : reagentStartWidths[level],
      isHorizontal ? reagentStartHeights[level] : levelStartY,
      level + 1,
      mol,
      isHorizontal,
    );

    if (isHorizontal)
      reagentStartHeights[level] += MOL_HEIGHT + VERTICAL_SPACING;
    else
      reagentStartWidths[level] += MOL_WIDTH + HORIZONTAL_SPACING;
  }
}

function reactionLinesToMoleculeSVG(
  svg: SVGElement,
  startX: number,
  startY: number,
  endX: number,
  endY: number,
  isHorizontal: boolean,
) {
  const group = document.createElementNS('http://www.w3.org/2000/svg', 'g');

  if (isHorizontal) {
    // Horizontal orientation: draw lines from right to left
    const startLineX = startX + MOL_WIDTH + PADDING / 2;
    const startLineY = startY + MOL_HEIGHT / 2;
    const endLineX = endX - PADDING / 2;
    const endLineY = endY + MOL_HEIGHT / 2;
    const midX = (startLineX + endLineX) / 2;

    drawLineSVG(group, startLineX, startLineY, midX, startLineY, true);
    drawLineSVG(group, midX, startLineY, midX, endLineY);
    drawLineSVG(group, midX, endLineY, endLineX, endLineY);
  } else {
    // Vertical orientation: draw lines from top to bottom
    const startLineX = startX + MOL_WIDTH / 2;
    const startLineY = startY + MOL_HEIGHT + PADDING / 2;
    const endLineX = endX + MOL_WIDTH / 2;
    const endLineY = endY - PADDING / 2;
    const midY = (startLineY + endLineY) / 2;

    drawLineSVG(group, startLineX, startLineY, startLineX, midY, true);
    drawLineSVG(group, startLineX, midY, endLineX, midY);
    drawLineSVG(group, endLineX, midY, endLineX, endLineY);
  }

  svg.appendChild(group);
}

async function drawMoleculeSVG(svg: SVGElement, smiles: string, x: number, y: number, color: string) {
  const group = document.createElementNS('http://www.w3.org/2000/svg', 'g');

  // Create clickable rectangle
  const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  setAttributes(rect, {'x': (x - PADDING / 2).toString(), 'y': (y - PADDING / 2).toString(),
    'width': (MOL_WIDTH + PADDING).toString(), 'height': (MOL_HEIGHT + PADDING).toString(), 'rx': '5',
    'stroke': color, 'fill': 'transparent', 'stroke-width': '1'});
  group.appendChild(rect);

  // Create foreignObject for the molecule rendering
  const foreignObject = document.createElementNS('http://www.w3.org/2000/svg', 'foreignObject');
  setAttributes(foreignObject, {'x': x.toString(), 'y': y.toString(), 'width': MOL_WIDTH.toString(),
    'height': MOL_HEIGHT.toString()});

  const molDiv = await grok.chem.drawMolecule(smiles, MOL_WIDTH, MOL_HEIGHT);
  molDiv.classList.add('retrosynthesis-mol-rect');

  foreignObject.appendChild(molDiv);
  group.appendChild(foreignObject);

  svg.appendChild(group);
}

function drawLineSVG(
  group: SVGGElement,
  startX: number,
  startY: number,
  endX: number,
  endY: number,
  withDot = false,
) {
  const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
  setAttributes(line, {'x1': startX.toString(), 'y1': startY.toString(), 'x2': endX.toString(), 'y2': endY.toString(),
    'stroke': 'black', 'stroke-width': '1'});
  group.appendChild(line);

  if (withDot) {
    const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    setAttributes(circle, {'cx': endX.toString(), 'cy': endY.toString(), 'r': '2', 'fill': 'black'});
    group.appendChild(circle);
  }
}

function setAttributes(element: SVGElement, attributes: {[key: string]: string}) {
  for (const key of Object.keys(attributes))
    element.setAttribute(key, attributes[key]);
}
