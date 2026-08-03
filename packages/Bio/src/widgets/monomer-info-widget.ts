/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';
import {MONOMER_MOTIF_SPLITTER, MONOMER_CANONICALIZER_TEMP} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {IMonomerCanonicalizer} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

import {getCorrectedSmiles, capSmiles} from '../utils/monomer-lib/monomer-manager/monomer-manager';

/** Caps the monomer (replaces R-groups with cap atoms) and returns capped SMILES. */
function getCappedSmiles(monomer: Monomer): string | null {
  try {
    const corrected = getCorrectedSmiles(monomer[REQ.RGROUPS], monomer.smiles, monomer.molfile);
    return capSmiles(corrected, monomer[REQ.RGROUPS]);
  } catch { /* capping may fail for some monomers */ }
  return null;
}

/** Renders a single-monomer details pane content. */
function renderMonomerDetails(monomer: Monomer): HTMLElement {
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
}

/** Renders the molecule info panel pane content for a single monomer. */
function renderMoleculePane(monomer: Monomer): HTMLElement {
  const cappedMol = getCappedSmiles(monomer);
  if (!cappedMol)
    return ui.divText('No molecular structure available');

  const molSv = DG.SemanticValue.fromValueType(cappedMol, DG.SEMTYPE.MOLECULE);
  return ui.panels.infoPanel(molSv).root;
}

/** Tries to get the canonicalizer from the SemanticValue's cell column. */
function getCanonicalizer(sv: DG.SemanticValue): IMonomerCanonicalizer | null {
  try {
    const col = sv.cell?.column;
    if (col)
      return col.temp[MONOMER_CANONICALIZER_TEMP] ?? null;
  } catch { /* no cell context */ }
  return null;
}

/** Creates a widget for the monomer info panel shown in the context panel. */
export function getMonomerInfoWidget(sv: DG.SemanticValue, monomerLib: IMonomerLib): DG.Widget {
  const rawValue = sv.value as string;
  if (!rawValue)
    return new DG.Widget(ui.divText('No monomer value.'));

  // Canonicalize if a canonicalizer is available from the column context
  const canonicalizer = getCanonicalizer(sv);
  const canonicalized = canonicalizer ? canonicalizer.canonicalize(rawValue) : rawValue;

  // Split by motif splitter — a cell value may contain multiple monomers
  const symbols = canonicalized.split(MONOMER_MOTIF_SPLITTER).map((s) => s.trim()).filter((s) => s.length > 0);

  // Resolve monomers from library
  const monomers: {symbol: string; monomer: Monomer | null}[] =
    symbols.map((s) => ({symbol: s, monomer: monomerLib.getMonomer(null, s)}));

  const found = monomers.filter((m) => m.monomer !== null);
  if (found.length === 0)
    return new DG.Widget(ui.divText(`Monomer '${rawValue}' not found in the library.`));

  const acc = ui.accordion('Monomer');

  if (found.length === 1) {
    // Single monomer — flat panes
    const m = found[0].monomer!;
    acc.addPane('Details', () => renderMonomerDetails(m), true);
    acc.addPane('Molecule', () => renderMoleculePane(m), true);
  } else {
    // Multiple monomers — one pane per monomer
    for (const {symbol, monomer} of found) {
      acc.addPane(symbol, () => {
        return ui.divV([
          renderMonomerDetails(monomer!),
          ui.element('hr'),
          renderMoleculePane(monomer!),
        ]);
      }, true);
    }
  }

  return new DG.Widget(acc.root);
}
