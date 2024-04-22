import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {printLeftOrCentered} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {findMonomers, parseHelm, removeGapsFromHelm} from './utils';
import {IEditor, HelmMonomerPlacer, IEditorMolAtom, ISeqMonomer, Temps} from './helm-monomer-placer';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';

// import {_package} from './package'; // NullError: method not found: '_package' on null

const enum tempTAGS {
  helmSumMaxLengthWords = 'helm-sum-maxLengthWords',
  helmMaxLengthWords = 'helm-maxLengthWords',

  helmPlacer = 'bio-helmPlacer',
}

function getSeqMonomerFromHelm(
  helmPrefix: string, symbol: string, monomerLib: IMonomerLib
): ISeqMonomer {
  let resSeqMonomer: ISeqMonomer | undefined = undefined;
  for (const polymerType of monomerLib.getPolymerTypes()) {
    if (helmPrefix.startsWith(polymerType))
      resSeqMonomer = {symbol: symbol, polymerType: polymerType};
  }
  if (!resSeqMonomer)
    throw new Error(`Monomer not found for symbol = '${symbol}' and helmPrefix = '${helmPrefix}'.`);
  return resSeqMonomer;
}

function getSeqMonomerFromHelmAtom(atom: IEditorMolAtom): ISeqMonomer {
  let polymerType: string | undefined = undefined;
  // @ts-ignore
  switch (atom.bio.type) {
  case 'HELM_BASE':
  case 'HELM_SUGAR': // r - ribose, d - deoxyribose
  case 'HELM_LINKER': // p - phosphate
    polymerType = 'RNA';
    break;
  case 'HELM_AA':
    polymerType = 'PEPTIDE';
    break;
  default:
    polymerType = 'PEPTIDE';
  }
  return {symbol: atom.elem, polymerType: polymerType};
}

function getHoveredMonomerFromEditor(
  argsX: number, argsY: number, gridCell: DG.GridCell, editor: IEditor
): ISeqMonomer | null {
  let hoveredSeqMonomer: ISeqMonomer | null = null;

  /** @return {[number, number]} [atom, distance] */
  function getNearest(excluded: (number | undefined)[]): [number | undefined, number | undefined] {
    let atom: number | undefined = undefined;
    let distance: number | undefined = undefined;
    for (let atomI = 0; atomI < editor.m.atoms.length; ++atomI) {
      if (!excluded.includes(atomI)) {
        const aX: number = editor.m.atoms[atomI].p.x;
        const aY: number = editor.m.atoms[atomI].p.y;
        const distanceToAtomI: number = Math.sqrt((argsX - aX) ** 2 + (argsY - aY) ** 2);
        if (distance === undefined || distance > distanceToAtomI) {
          atom = atomI;
          distance = distanceToAtomI;
        }
      }
    }
    return [atom, distance];
  }

  const [firstAtomI, firstDistance] = getNearest([]);
  const [secondAtomI, secondDistance] = getNearest([firstAtomI]);

  // for (let atomI = 0; atomI < editor.m.atoms.length; ++atomI) {
  //   const aX: number = editor.m.atoms[atomI].p.x;
  //   const aY: number = editor.m.atoms[atomI].p.y;
  //   const distanceToAtomI: number = Math.sqrt((argsX - aX) ** 2 + (argsY - aY) ** 2);
  //   if (firstDistance === undefined || firstDistance > distanceToAtomI) {
  //     secondAtomI = firstAtomI;
  //     firstAtomI = atomI;
  //
  //     secondDistance = firstDistance;
  //     firstDistance = distanceToAtomI;
  //   }
  // }

  if (firstAtomI !== undefined && firstDistance !== undefined) {
    const firstAtom = editor.m.atoms[firstAtomI];
    const firstSeqMonomer = getSeqMonomerFromHelmAtom(firstAtom);
    if (secondAtomI !== undefined && secondDistance !== undefined) {
      if (firstDistance < secondDistance * 0.45)
        hoveredSeqMonomer = firstSeqMonomer;
    } else {
      if (firstDistance < 0.35 * gridCell.bounds.height)
        hoveredSeqMonomer = firstSeqMonomer;
    }
  }
  return hoveredSeqMonomer;
}

function getHoveredMonomerFallback(
  argsX: number, _argsY: number, gridCell: DG.GridCell, helmPlacer: HelmMonomerPlacer
): ISeqMonomer | null {
  let hoveredSeqMonomer: ISeqMonomer | null = null;
  const [allParts, lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);
  const maxIndex = Object.values(lengths).length - 1;
  let left = 0;
  let right = maxIndex;
  let found = false;
  let iterCount: number = 0;

  let mid = 0;
  if (argsX > sumLengths[0]) {
    while (!found && iterCount < sumLengths.length) {
      mid = Math.floor((right + left) / 2);
      if (argsX >= sumLengths[mid] && argsX <= sumLengths[mid + 1]) {
        left = mid;
        found = true;
      } else if (argsX < sumLengths[mid])
        right = mid - 1;
      else if (argsX > sumLengths[mid + 1])
        left = mid + 1;

      if (left == right)
        found = true;

      iterCount++;
    }
  }
  left = (argsX >= sumLengths[left]) ? left : left - 1; // correct left to between sumLengths
  if (left >= 0)
    hoveredSeqMonomer = getSeqMonomerFromHelm(allParts[0], allParts[left], helmPlacer.monomerLib);
  return hoveredSeqMonomer;
}

type RendererGridCellTemp = {
  [Temps.helmMonomerPlacer]: HelmMonomerPlacer;
}

function getRendererGridCellTemp(gridCell: DG.GridCell
): [DG.GridColumn | null, DG.Column | null, RendererGridCellTemp] {
  /** Primarily store/get MonomerPlacer at GridColumn, fallback at (Table) Column for scatter plot tooltip  */
  let temp: RendererGridCellTemp | null = null;

  let gridCol: DG.GridColumn | null = null;
  try { gridCol = gridCell.gridColumn; } catch { gridCol = null; }
  temp = gridCol ? gridCol.temp as RendererGridCellTemp : null;

  const tableCol: DG.Column = gridCell.cell.column;
  if (!temp) temp = tableCol.temp;

  if (temp === null) throw new Error(`Monomer placer store (GridColumn or Column) not found.`);
  return [gridCol, tableCol, temp];
}

/** Helm cell renderer in case of no missed monomer draws with JSDraw2.Editor (webeditor),
 * in case of missed monomers presented, draws linear sequences aligned in width per monomer.
 */
export class HelmCellRenderer extends DG.GridCellRenderer {
  get name() { return 'helm'; }

  get cellType() { return 'helm'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    try {
      const [_gridCol, tableCol, temp] = getRendererGridCellTemp(gridCell);
      const helmPlacer = temp[Temps.helmMonomerPlacer];
      /* Can not do anything without tableColumn */
      if (!tableCol) return;

      /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
      const argsX = e.offsetX - gcb.x;
      const argsY = e.offsetY - gcb.y;

      const editor: IEditor | null = helmPlacer.getEditor(gridCell.tableRowIndex!);
      let seqMonomer: ISeqMonomer | null;
      let missedMonomers: Set<string> = new Set<string>(); // of .size = 0
      if (editor)
        seqMonomer = getHoveredMonomerFromEditor(argsX, argsY, gridCell, editor);
      else {
        const seq: string = !gridCell.cell.value ? '' : removeGapsFromHelm(gridCell.cell.value as string);
        const monomerList = parseHelm(seq);
        missedMonomers = findMonomers(monomerList);
        const parsedMonomers = new Set<string>(monomerList);
        seqMonomer = getHoveredMonomerFallback(argsX, argsY, gridCell, helmPlacer);
        if (seqMonomer && !parsedMonomers.has(seqMonomer.symbol)) seqMonomer = null;
        if (seqMonomer) {
          const textSize = helmPlacer.monomerTextSizeMap[seqMonomer.symbol];
          if (textSize) {
            const textBaseLine = gcb.height / 2 -
              (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2 + 1;
            const textTop = textBaseLine - textSize.fontBoundingBoxAscent;
            const textBottom = textBaseLine + textSize.fontBoundingBoxDescent;
            if (argsY < textTop || textBottom < argsY) seqMonomer = null;
          }
        }
      }

      if (seqMonomer) {
        if (!missedMonomers.has(seqMonomer.symbol)) {
          const tooltipElements: HTMLElement[] = [ui.div(seqMonomer.symbol)];
          const monomerDiv = helmPlacer.monomerLib.getTooltip(seqMonomer.polymerType, seqMonomer.symbol);
          tooltipElements.push(monomerDiv);
          ui.tooltip.show(ui.divV(tooltipElements), e.x + 16, e.y + 16);
        } else {
          ui.tooltip.show(ui.divV([
            ui.divText(`Monomer '${seqMonomer.symbol}' not found.`),
            ui.divText('Open the Context Panel, then expand Manage Libraries'),
          ]), e.x + 16, e.y + 16);
        }
      } else if (missedMonomers.size == 0) {
        ui.tooltip.hide();
        return;
      } else { // seqMonomer == null && missedMonomers.size > 0
        const mmStrList = wu(missedMonomers.keys()).toArray().sort()
          .filter((_, i) => i < 3);
        const missedMonomersStr = mmStrList.join(', ') + (missedMonomers.size > 3 ? ', ...' : '');
        ui.tooltip.show(ui.divV([ui.divText('Monomers missed in monomer libraries:'),
          ui.divText(missedMonomersStr)
        ]), e.x + 16, e.y + 16);
      }
    } catch (err: any) {
      const errMsg: string = errorToConsole(err);
      console.error('Helm: HelmCellRenderer.onMouseMove() error:\n' + errMsg);
    } finally {
      e.preventDefault();
      e.stopPropagation();
    }
  }


  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    g.save();
    try {
      const [gridCol, tableCol, temp] = getRendererGridCellTemp(gridCell);
      /* Can not do anything without tableColumn containing temp */
      if (!tableCol) return;
      let helmPlacer = temp[Temps.helmMonomerPlacer];
      if (!helmPlacer) helmPlacer = temp[Temps.helmMonomerPlacer] = new HelmMonomerPlacer(gridCol, tableCol);

      const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
      const missedColor = 'red';
      const monomerColor: string = '#404040';
      const frameColor: string = '#C0C0C0';

      const seq: string = !gridCell.cell.value ? '' : removeGapsFromHelm(gridCell.cell.value);
      const monomerList = parseHelm(seq);
      const monomers: Set<string> = new Set<string>(monomerList);
      const missedMonomers: Set<string> = findMonomers(monomerList);

      if (missedMonomers.size == 0) {
        // Recreate host to avoid hanging in window.dojox.gfx.svg.Text.prototype.getTextWidth
        const host = gridCell.element = ui.div([],
          {style: {width: `${w - 2}px`, height: `${h - 2}px`, margin: `${1}px`, backgroundColor: '#FFE0E0'}});
        host.setAttribute('dataformat', 'helm');
        host.setAttribute('data', seq /* gaps skipped */);
        // if grid has neighbour to the left, then shift host to the left
        if (host.parentElement && (gridCell.grid?.canvas?.offsetLeft ?? 0) > 0) {
          host.parentElement.style.left =
            `${(gridCell.grid?.canvas?.offsetLeft ?? 0) + host.parentElement.offsetLeft}px`;
        }

        // Recreate editor to avoid hanging in window.dojox.gfx.svg.Text.prototype.getTextWidth
        const editor = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true}) as IEditor;
        helmPlacer.setEditor(gridCell.tableRowIndex!, editor);

        helmPlacer.skipCell(gridCell.tableRowIndex!);
        return;
      }

      if (missedMonomers.size > 0) {
        if (!grid) {
          const r = window.devicePixelRatio;
          h = 28;
          g.canvas.height = h * r;
          g.canvas.style.height = `${h}px`;
        }

        w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
        //g.save();
        g.beginPath();
        g.rect(x, y, w, h);
        g.clip();
        g.transform(1, 0, 0, 1, x, y);
        g.font = '12px monospace';
        g.textBaseline = 'top';
        const [allParts, _lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);

        for (let i = 0; i < allParts.length; ++i) {
          const part: string = allParts[i];
          const color: string =
            part === '.' || part.endsWith('{') || part.startsWith('}') ? frameColor :
              missedMonomers.has(part) ? missedColor :
                monomers.has(part) ? monomerColor :
                  frameColor;
          g.fillStyle = color;
          printLeftOrCentered(sumLengths[i], 0, w, h, g, allParts[i], color, 0, true, 1.0,
            undefined, undefined, undefined, undefined, undefined,
            undefined, undefined, undefined, helmPlacer.monomerTextSizeMap);
        }
      }
    } finally {
      g.restore();
    }
  }
}
