import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {byEngineModule} from '../engine/chem-engine';
import {noteCruxError} from '../engine/crux-module';

const WIDTH = 200;
const HEIGHT = 100;
const MORE_ICON_FONT_WEIGHT = '500';

export function structure2dWidget(molecule: string): DG.Widget {
  const readable = byEngineModule('structure2d', getRdKitModule(), (module, crux) => {
    const mol = getMolSafe(molecule, {}, module, true, crux ? noteCruxError : undefined).mol;
    mol?.delete();
    return mol ? true : crux ? null : false;
  });
  const resultWidget = readable ?
    new DG.Widget(get2dMolecule(molecule, WIDTH, HEIGHT)) : new DG.Widget(ui.divText('Molecule is possibly malformed'));
  //in case we undock the widget - redraw the molecule when resize
  ui.tools.waitForElementInDom(resultWidget.root).then(() => {
    if (resultWidget.root.closest('.dialog-floating')) {
      const accPanel = resultWidget.root.closest('.panel-content') as HTMLElement;
      if (accPanel) {
        ui.onSizeChanged(accPanel).subscribe((_) => {
          const w = Math.max(accPanel.clientWidth, 200);
          const h = Math.max(accPanel.clientHeight, 100);
          const waitParentEl = resultWidget.root.parentElement;
          if (waitParentEl?.classList.contains('grok-wait')) {
            waitParentEl.style.width = `${w}px`;
            waitParentEl.style.height = `${h}px`;
          }
          ui.empty(resultWidget.root);
          resultWidget.root.append(get2dMolecule(molecule, w, h, true));
        });
      }
    }
  });
  return resultWidget;
}

function get2dMolecule(molecule: string, w: number, h: number, undocked?: boolean): HTMLElement {
  const molecule2d = renderMolecule(molecule, {renderer: 'RDKit', width: w, height: h});
  const moreIcon = molecule2d.querySelector('.chem-mol-view-icon.pep-more-icon') as HTMLElement;
  undocked ? moreIcon.style.right = `50px` : moreIcon.style.left = `${w}px`;
  moreIcon.style.fontWeight = MORE_ICON_FONT_WEIGHT;
  moreIcon.classList.remove('pep-more-icon');
  molecule2d.classList.remove('d4-flex-col');
  molecule2d.classList.add('d4-flex-wrap');
  return ui.div(molecule2d, 'd4-flex-col');
}
