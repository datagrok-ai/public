import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parsePdbHeaders} from './pdb-parser';

/** Fetch SMILES for a list of ligand compound IDs from RCSB Chemical Component Dictionary. */
async function fetchLigandSmiles(compIds: string[]): Promise<{[compId: string]: string}> {
  const query = `{
    chem_comps(comp_ids: [${compIds.map((id) => `"${id}"`).join(', ')}]) {
      rcsb_id
      rcsb_chem_comp_descriptor {
        SMILES_stereo
        SMILES
      }
    }
  }`;
  try {
    const resp = await grok.dapi.fetchProxy('https://data.rcsb.org/graphql', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({query}),
    });
    if (!resp.ok) return {};
    const data = await resp.json();
    const result: {[compId: string]: string} = {};
    for (const comp of data.data?.chem_comps ?? []) {
      const smiles = comp.rcsb_chem_comp_descriptor?.SMILES_stereo ??
        comp.rcsb_chem_comp_descriptor?.SMILES;
      if (smiles)
        result[comp.rcsb_id] = smiles;
    }
    return result;
  } catch {
    return {};
  }
}

export function pdbFileInfoWidget(pdbText: string): DG.Widget {
  const info = parsePdbHeaders(pdbText);

  // -- Basic info table --
  const basicMap: {[key: string]: any} = {};
  if (info.title)
    basicMap['Description'] = info.title;
  if (info.classification)
    basicMap['Classification'] = info.classification;
  if (info.method)
    basicMap['Method'] = info.method;
  if (info.resolution != null)
    basicMap['Resolution'] = `${info.resolution} \u00C5`;
  if (info.pdbId) {
    basicMap['PDB URL'] = ui.link(
      `rcsb.org/structure/${info.pdbId}`,
      () => window.open(`https://www.rcsb.org/structure/${info.pdbId}`),
    );
  }
  if (info.organism)
    basicMap['Organism'] = info.organism;
  if (info.expressionSystem)
    basicMap['Expression System'] = info.expressionSystem;
  if (info.atomCount != null)
    basicMap['Atom Count'] = info.atomCount.toString();
  if (info.chains)
    basicMap['Chains'] = info.chains.join(', ');

  // If absolutely nothing was parsed, show a minimal message
  if (Object.keys(basicMap).length === 0)
    return new DG.Widget(ui.divText('No header information available'));

  const basicTable = ui.tableFromMap(basicMap);

  // -- Accordion for subsections --
  const acc = DG.Accordion.create();

  // -- Structure Details --
  const hasDetails = info.rFactor != null || info.disulfideBondCount != null || info.depositDate;
  if (hasDetails) {
    acc.addPane('Structure Details', () => {
      const map: {[key: string]: any} = {};
      if (info.depositDate)
        map['Deposit Date'] = info.depositDate;
      if (info.rFactor != null)
        map['R-factor'] = info.rFactor.toFixed(3);
      if (info.disulfideBondCount != null)
        map['Disulfide Bonds'] = info.disulfideBondCount.toString();
      return ui.tableFromMap(map);
    }, false);
  }

  // -- Literature --
  if (info.citation) {
    acc.addPane('Literature', () => {
      const citMap: {[key: string]: any} = {};
      if (info.citation!.title)
        citMap['Name'] = info.citation!.title;
      if (info.citation!.journal)
        citMap['Journal'] = info.citation!.journal;
      if (info.citation!.doi) {
        citMap['DOI'] = ui.link(
          info.citation!.doi,
          () => window.open(`https://doi.org/${info.citation!.doi}`),
        );
      }
      if (info.citation!.pmid) {
        citMap['PubMed'] = ui.link(
          info.citation!.pmid.toString(),
          () => window.open(`https://pubmed.ncbi.nlm.nih.gov/${info.citation!.pmid}`),
        );
      }
      return ui.tableFromMap(citMap);
    }, false);
  }

  // -- Macromolecules --
  if (info.entities && info.entities.length > 0) {
    acc.addCountPane('Macromolecules', () => {
      const innerAcc = DG.Accordion.create();
      for (const entity of info.entities!) {
        const label = `Entity ${entity.molId}` + (entity.molecule ? `: ${entity.molecule}` : '');
        innerAcc.addPane(label, () => {
          const map: {[key: string]: any} = {};
          if (entity.chains)
            map['Chains'] = entity.chains.join(', ');
          if (entity.molecule)
            map['Molecule'] = entity.molecule;
          return ui.tableFromMap(map);
        }, false);
      }
      return innerAcc.root;
    }, () => info.entities!.length, false);
  }

  // -- Small Molecules (Ligands) --
  if (info.ligands && info.ligands.length > 0) {
    // Determine if this is an RCSB PDB (has valid 4-char PDB ID) for SMILES fetching
    const isRcsb = !!(info.pdbId && /^[A-Za-z0-9]{4}$/.test(info.pdbId));

    acc.addCountPane('Small Molecules', () => {
      const innerAcc = DG.Accordion.create();
      for (const lig of info.ligands!) {
        innerAcc.addPane(lig.id, () => {
          const map: {[key: string]: any} = {};
          if (lig.chain)
            map['Chain'] = lig.chain;
          if (lig.name)
            map['Name'] = lig.name;
          if (lig.formula)
            map['Formula'] = lig.formula;
          const parts: HTMLElement[] = [ui.tableFromMap(map)];

          // Fetch 2D structure from RCSB if this is a standard PDB
          if (isRcsb) {
            const molHost = ui.div([ui.loader()]);
            parts.push(molHost);
            fetchLigandSmiles([lig.id]).then((smilesMap) => {
              molHost.innerHTML = '';
              const smiles = smilesMap[lig.id];
              if (smiles)
                molHost.appendChild(grok.chem.drawMolecule(smiles, 250, 200));
            });
          }

          return ui.divV(parts);
        }, false);
      }
      return innerAcc.root;
    }, () => info.ligands!.length, false);
  }

  const parts: HTMLElement[] = [basicTable];
  if (acc.panes.length > 0)
    parts.push(acc.root);

  return new DG.Widget(ui.divV(parts));
}
