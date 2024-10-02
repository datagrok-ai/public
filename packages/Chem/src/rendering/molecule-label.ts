import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../package';
import {RDKIT_COMMON_RENDER_OPTS} from '../utils/chem-common-rdkit';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getFirstNSymbols} from '../utils/chem-common';

export async function drawMoleculeLabels(
  table: DG.DataFrame, molCol: DG.Column, sp: DG.ScatterPlotViewer, maxPoints: number, smallMarkerSize: number = -1,
  largeMarkerSize: number, minDistBetweenLabels: number, sizeCol: string = ''): Promise<void> {
  let smallMarker = true;

  let smallMarkerType = sp.getOptions().look['markerType'];
  sp.onAfterDrawScene.pipe().subscribe(() => {
    if (smallMarker)
      smallMarkerType = sp.getOptions().look['markerType'];
    const rowCount = table.rowCount;
    const pointsOnScreen = new Array<DG.Point>(maxPoints);
    const pointsOnScreenIdxs = new Uint32Array(maxPoints);
    const {maxX, maxY, minX, minY} = sp.viewBox;

    let counter = 0;
    for (let i = 0; i < rowCount; i++) {
      //@ts-ignore
      const point = sp.pointToScreen(i);
      if (point.x > minX && point.x < maxX && point.y > minY && point.y < maxY && table.filter.get(i)) {
        if (counter == 20) {
          sp.setOptions({
            markerDefaultSize: smallMarkerSize,
            sizeColumnName: sizeCol,
            markerType: smallMarkerType,
          });
          if (!smallMarker) {
            smallMarker = true;
            setTimeout(() => sp.render(sp.getInfo()['canvas'].getContext('2d')), 10);
          }
          return; // return if there are more than max allowed points in total on the screen
        }
        else {
          pointsOnScreen[counter] = point;
          pointsOnScreenIdxs[counter] = i;
          counter++;
        }
      }
    }

    const N = counter * (counter + 1) / 2;
    const distancesMatrix = new Float32Array(N);

    for (let i = 0; i < counter; i++) {
      for (let j = 0; j < counter; j++) {
        if (i === j)
          continue;
        else {
          const dist = j < i ? distancesMatrix[(j * counter - (j - 1) * j / 2) + (i - j) - 1] :
            pointsOnScreen[i]!.distanceTo(pointsOnScreen[j]!);
          if (j > i)
            distancesMatrix[(i * counter - (i - 1) * i / 2) + (j - i) - 1] = dist;
          if (dist < minDistBetweenLabels) {
            sp.setOptions({
              markerDefaultSize: smallMarkerSize,
              sizeColumnName: sizeCol,
              markerType: smallMarkerType,
            });
            if (!smallMarker) {
              smallMarker = true;
              setTimeout(() => sp.render(sp.getInfo()['canvas'].getContext('2d')), 10);
            }
            return; // return if distance between points is less than min allowed
          }
        }
      }
    }

    const ctxMain = sp.getInfo()['canvas'].getContext('2d') as CanvasRenderingContext2D;
    //const sizeColName = sizeCol ? sizeCol.startsWith('sali') ? 'sali' : sizeCol : '';
    ctxMain.textAlign = 'left';
    //ctxMain.textBaseline = 'top';
    const textOffset = 20;
    const rdkitModule = getRdKitModule();

    for (let i = 0; i < counter; i++) {
      let mol: RDMol | null = null;
      try {
        ctxMain.beginPath();
        ctxMain.strokeStyle = `rgba(0, 0, 0, 0.1)`;
        ctxMain.fillStyle = `white`;
        const moleculeRectOffset = 10;
        const canwasWidth = largeMarkerSize - moleculeRectOffset;
        const canvasHeight = canwasWidth * 0.7;
        const imageHost = ui.canvas(canwasWidth, canvasHeight);
        mol = rdkitModule.get_mol(molCol.get(pointsOnScreenIdxs[i]));
        mol.draw_to_canvas_with_highlights(imageHost, JSON.stringify(RDKIT_COMMON_RENDER_OPTS));
        ctxMain.arc(pointsOnScreen[i].x, pointsOnScreen[i].y, Math.floor(largeMarkerSize/2), 0, 2 * Math.PI);
        ctxMain.fill();
        const leftX = pointsOnScreen[i].x - Math.floor(canwasWidth / 2);
        ctxMain.drawImage(imageHost, leftX, pointsOnScreen[i].y - Math.floor(canvasHeight / 2));
        if (sizeCol) {
          ctxMain.fillStyle = `black`;
          ctxMain.fillText(`${sizeCol}: ${getFirstNSymbols(table.get(sizeCol, pointsOnScreenIdxs[i]), 4)}`,
            leftX + textOffset, pointsOnScreen[i].y + Math.floor(canvasHeight / 2));
        }
        ctxMain.closePath();
      } finally {
        mol?.delete();
      }
    }
    sp.setOptions({
      sizeColumnName: '',
      markerDefaultSize: largeMarkerSize,
      markerType: 'circle border',
    });
    if (smallMarker) {
      smallMarker = false;
      setTimeout(() => sp.render(sp.getInfo()['canvas'].getContext('2d')), 10);
    }
  });
}
