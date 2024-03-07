import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {printLeftOrCentered} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {findMonomers, parseHelm} from './utils';
import {IEditor, HelmMonomerPlacer, IEditorMolAtom, ISeqMonomer} from './helm-monomer-placer';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';

// import {_package} from './package'; // NullError: method not found: '_package' on null

const enum tempTAGS {
  helmSumMaxLengthWords = 'helm-sum-maxLengthWords',
  helmMaxLengthWords = 'helm-maxLengthWords',

  helmPlacer = 'bio-helmPlacer',
}

// Global flag is for replaceAll
const helmGapStartRe = /\{(\*\.)+/g;
const helmGapIntRe = /\.(\*\.)+/g;
const helmGapEndRe = /(\.\*)+\}/g;

function getSeqMonomerFromHelm(
  helmPrefix: string, symbol: string, monomerLib: IMonomerLib
): ISeqMonomer {
  let resSeqMonomer: ISeqMonomer | undefined = undefined;
  for (const polymerType of monomerLib.getPolymerTypes()) {
    if (helmPrefix.startsWith(polymerType))
      resSeqMonomer = {symbol: symbol, polymerType: polymerType};
  }
  if (!resSeqMonomer)
    resSeqMonomer = {symbol: symbol};
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
      /* Can not do anything without tableColumn containing temp */
      let tableCol: DG.Column | null = null;
      try { tableCol = gridCell.tableColumn; } catch { }
      if (!tableCol) return;

      const argsX = e.offsetX - gridCell.bounds.x;
      const argsY = e.offsetY - gridCell.bounds.y;

      const helmPlacer = HelmMonomerPlacer.getOrCreate(tableCol);
      const editor: IEditor | null = helmPlacer.getEditor(gridCell.tableRowIndex!);
      const seqMonomer: ISeqMonomer | null = editor ? getHoveredMonomerFromEditor(argsX, argsY, gridCell, editor) :
        getHoveredMonomerFallback(argsX, argsY, gridCell, helmPlacer);
      if (!seqMonomer) {
        ui.tooltip.hide();
        return;
      }

      const seq: string = !gridCell.cell.value ? '' : gridCell.cell.value
        .replaceAll(helmGapStartRe, '{').replaceAll(helmGapIntRe, '.').replaceAll(helmGapEndRe, '}')
        .replace('{*}', '{}');
      const monomerList = parseHelm(seq);
      const missedMonomers = findMonomers(monomerList);

      if (missedMonomers.has(seqMonomer.symbol)) {
        ui.tooltip.show(ui.divV([
          ui.divText(`Monomer ${seqMonomer.symbol} not found.`),
          ui.divText('Open the Context Panel, then expand Manage Libraries'),
        ]), e.x + 16, e.y + 16);
      } else {
        const monomer = helmPlacer.getMonomer(seqMonomer);
        if (monomer) {
          const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
          const monomerSvg = grok.chem.svgMol(monomer.smiles, undefined, undefined, options);
          ui.tooltip.show(ui.divV([
            ui.divText(seqMonomer.symbol),
            monomerSvg,
          ]), e.x + 16, e.y + 16);
        }
      }

      // const tooltipMessage: HTMLElement[] = [];
      // for (const [part, partI] of wu.enumerate(allParts)) {
      //   if (missedMonomers.has(part)) {
      //     tooltipMessage[partI] = ui.divV([
      //       ui.divText(`Monomer ${allParts[partI]} not found.`),
      //       ui.divText('Open the Context Panel, then expand Manage Libraries')
      //     ]);
      //   } else if (monomers.has(part)) {
      //     const elList = [ui.div(part)];
      //     const monomer = helmPlacer.getMonomer(part);
      //     if (monomer) {
      //       const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      //       const monomerSvg = grok.chem.svgMol(monomer.smiles, undefined, undefined, options);
      //       elList.push(monomerSvg);
      //     }
      //     tooltipMessage[partI] = ui.divV(elList);
      //   }
      // }
      //
      // (((tooltipMessage[left]?.childNodes.length ?? 0) > 0)) ?
      //   ui.tooltip.show(tooltipMessage[left], e.x + 16, e.y + 16) :
      //   ui.tooltip.hide();
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
      /* Can not do anything without tableColumn containing temp */
      let tableCol: DG.Column | null = null;
      try { tableCol = gridCell.tableColumn; } catch { }
      if (!tableCol) return;

      const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
      const missedColor = 'red';
      const monomerColor: string = '#404040';
      const frameColor: string = '#C0C0C0';

      const seq: string = !gridCell.cell.value ? '' : gridCell.cell.value
        .replaceAll(helmGapStartRe, '{').replaceAll(helmGapIntRe, '.').replaceAll(helmGapEndRe, '}')
        .replace('{*}', '{}');
      const monomerList = parseHelm(seq);
      const monomers: Set<string> = new Set<string>(monomerList);
      const missedMonomers: Set<string> = findMonomers(monomerList);
      const helmPlacer = HelmMonomerPlacer.getOrCreate(tableCol);

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
        const [allParts, lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);

        for (let i = 0; i < allParts.length; ++i) {
          const part: string = allParts[i];
          const color: string =
            part === '.' || part.endsWith('{') || part.startsWith('}') ? frameColor :
              missedMonomers.has(part) ? missedColor :
                monomers.has(part) ? monomerColor :
                  frameColor;
          g.fillStyle = color;
          printLeftOrCentered(sumLengths[i], 0, w, h, g, allParts[i], color, 0, true, 1.0);
        }
      }
    } finally {
      g.restore();
    }
  }
}
