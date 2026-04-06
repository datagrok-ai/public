/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';

import {getCorrectedSmiles, capSmiles} from '../utils/monomer-lib/monomer-manager/monomer-manager';

/** Finds a monomer in the library across all polymer types. */
function findMonomer(monomerLib: IMonomerLib, symbol: string): Monomer | null {
  return monomerLib.getMonomer(null, symbol);
}

/** Caps the monomer (replaces R-groups with cap atoms) and returns capped SMILES. */
function getCappedSmiles(monomer: Monomer): string | null {
  try {
    const corrected = getCorrectedSmiles(monomer[REQ.RGROUPS], monomer.smiles, monomer.molfile);
    return capSmiles(corrected, monomer[REQ.RGROUPS]);
  } catch { /* capping may fail for some monomers */ }
  return null;
}

/** Creates a widget for the monomer info panel shown in the context panel. */
export function getMonomerInfoWidget(symbol: string, monomerLib: IMonomerLib): DG.Widget {
  const monomer = findMonomer(monomerLib, symbol);
  if (!monomer)
    return new DG.Widget(ui.divText(`Monomer '${symbol}' not found in the library.`));

  const acc = ui.accordion('Monomer');

  // Details pane — includes general info, structure, and R-groups in one table
  acc.addPane('Details', () => {
    const map: {[key: string]: any} = {
      'Symbol': monomer[REQ.SYMBOL],
      'Name': monomer[REQ.NAME],
      'Polymer Type': monomer[REQ.POLYMER_TYPE],
      'Monomer Type': monomer[REQ.MONOMER_TYPE],
    };
    if (monomer[REQ.AUTHOR])
      map['Author'] = monomer[REQ.AUTHOR];
    if (monomer.naturalAnalog)
      map['Natural Analog'] = monomer.naturalAnalog;
    if (monomer.lib?.source) {
      let source = monomer.lib.source;
      if (source.endsWith('.json'))
        source = source.substring(0, source.length - 5);
      map['Library'] = source;
    }

    // Structure
    if (monomer.molfile)
      map['Structure'] = grok.chem.drawMolecule(monomer.molfile, 150, 150);
    else if (monomer.smiles)
      map['Structure'] = grok.chem.drawMolecule(monomer.smiles, 150, 150);

    // R-Groups
    const rgroups = monomer[REQ.RGROUPS];
    if (rgroups && rgroups.length > 0) {
      for (const rg of rgroups)
        map[rg[RGP.LABEL] ?? '?'] = rg[RGP.ALTERNATE_ID] ?? '';
    }

    return ui.tableFromMap(map);
  }, true);

  // Molecule panel pane — cap the monomer first, then embed the generic Molecule context panel
  acc.addPane('Molecule', () => {
    const cappedMol = getCappedSmiles(monomer);
    if (!cappedMol)
      return ui.divText('No molecular structure available');

    const molSv = DG.SemanticValue.fromValueType(cappedMol, DG.SEMTYPE.MOLECULE);
    return ui.panels.infoPanel(molSv).root;
  }, true);

  return new DG.Widget(acc.root);
}
