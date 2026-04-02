import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parsePdbHeaders} from './pdb-helper';

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
  if (info.softwareRemarks)
    basicMap['Prepared By'] = info.softwareRemarks.join('; ');
  if (info.atomCount != null)
    basicMap['Atom Count'] = info.atomCount.toString();
  if (info.chains)
    basicMap['Chains'] = info.chains.join(', ');

  // If absolutely nothing was parsed, show a minimal message
  if (Object.keys(basicMap).length === 0)
    return new DG.Widget(ui.divText('No header information available'));

  const basicTable = ui.divV([ui.h2('General'), ui.tableFromMap(basicMap)]);

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
          // Show residue count per chain
          if (entity.chains && info.residueCountByChain) {
            const counts = entity.chains
              .filter((c) => info.residueCountByChain![c] != null)
              .map((c) => `${c}: ${info.residueCountByChain![c]}`);
            if (counts.length > 0)
              map['Residues'] = counts.join(', ');
          }
          if (entity.engineered)
            map['Engineered'] = 'Yes';
          if (entity.mutation)
            map['Mutation'] = 'Yes';
          // Show specific mutations for this entity's chains
          if (info.mutations && entity.chains) {
            const entityMuts = info.mutations.filter((m) => entity.chains!.includes(m.chain));
            if (entityMuts.length > 0)
              map['Mutation(s)'] = entityMuts.map((m) => m.description).join(', ');
          }
          return ui.tableFromMap(map);
        }, false);
      }
      return innerAcc.root;
    }, () => info.entities!.length, false);
  }

  // Determine if this is an RCSB PDB (has valid 4-char PDB ID) for SMILES fetching
  const isRcsb = !!(info.pdbId && /^[A-Za-z0-9]{4}$/.test(info.pdbId));

  // -- Small Molecules (Ligands) --
  if (info.ligands && info.ligands.length > 0) {
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
          const tableHost = ui.div();
          const parts: HTMLElement[] = [tableHost];

          // Fetch SMILES and 2D structure from RCSB if this is a standard PDB
          if (isRcsb) {
            const molHost = ui.div([ui.loader()]);
            parts.push(molHost);
            tableHost.appendChild(ui.tableFromMap(map));
            fetchLigandSmiles([lig.id]).then((smilesMap) => {
              const smiles = smilesMap[lig.id];
              if (smiles) {
                map['SMILES'] = smiles;
                tableHost.innerHTML = '';
                tableHost.appendChild(ui.tableFromMap(map));
                molHost.innerHTML = '';
                molHost.appendChild(grok.chem.drawMolecule(smiles, 250, 200));
              } else
                molHost.innerHTML = '';
            });
          } else
            tableHost.appendChild(ui.tableFromMap(map));


          return ui.divV(parts);
        }, false);
      }
      return innerAcc.root;
    }, () => info.ligands!.length, false);
  }

  // -- Binding Sites --
  if (info.sites && info.sites.length > 0) {
    acc.addCountPane('Binding Sites', () => {
      const innerAcc = DG.Accordion.create();
      for (const site of info.sites!) {
        innerAcc.addPane(site.name, () => {
          return ui.divText(site.residues.join(', '));
        }, false);
      }
      return innerAcc.root;
    }, () => info.sites!.length, false);
  }

  // -- Secondary Structure --
  if (info.secondaryStructure) {
    const ss = info.secondaryStructure;
    acc.addPane('Secondary Structure', () => {
      const map: {[key: string]: any} = {};
      if (ss.helices.length > 0) {
        map['Helices'] = ss.helices.length.toString();
        // Summarize by chain
        const helixByChain: {[c: string]: number} = {};
        for (const h of ss.helices)
          helixByChain[h.chain] = (helixByChain[h.chain] || 0) + 1;
        map['Helices by Chain'] = Object.entries(helixByChain)
          .map(([c, n]) => `${c}: ${n}`).join(', ');
      }
      if (ss.sheets.length > 0) {
        map['Strands'] = ss.sheets.length.toString();
        const sheetByChain: {[c: string]: number} = {};
        for (const s of ss.sheets)
          sheetByChain[s.chain] = (sheetByChain[s.chain] || 0) + 1;
        map['Strands by Chain'] = Object.entries(sheetByChain)
          .map(([c, n]) => `${c}: ${n}`).join(', ');
      }
      return ui.tableFromMap(map);
    }, false);
  }

  // -- Modified Residues --
  if (info.modifiedResidues && info.modifiedResidues.length > 0) {
    acc.addCountPane('Modified Residues', () => {
      const innerAcc = DG.Accordion.create();
      for (const mod of info.modifiedResidues!) {
        const label = `${mod.resName} ${mod.chain}${mod.resSeq}`;
        innerAcc.addPane(label, () => {
          const map: {[key: string]: any} = {};
          map['Modified Residue'] = mod.resName;
          map['Standard Residue'] = mod.standardRes;
          map['Chain'] = mod.chain;
          map['Position'] = mod.resSeq.toString();
          if (mod.comment)
            map['Description'] = mod.comment;
          const modParts: HTMLElement[] = [ui.tableFromMap(map)];

          // Fetch 2D structure from RCSB if this is a standard PDB
          if (isRcsb) {
            const molHost = ui.div([ui.loader()]);
            modParts.push(molHost);
            fetchLigandSmiles([mod.resName]).then((smilesMap) => {
              molHost.innerHTML = '';
              const smiles = smilesMap[mod.resName];
              if (smiles)
                molHost.appendChild(grok.chem.drawMolecule(smiles, 250, 200));
            });
          }

          return ui.divV(modParts);
        }, false);
      }
      return innerAcc.root;
    }, () => info.modifiedResidues!.length, false);
  }

  // -- Missing Residues --
  if (info.missingResidues && info.missingResidues.length > 0) {
    acc.addCountPane('Missing Residues', () => {
      // Group by chain
      const byChain: {[chain: string]: Array<{resName: string; resSeq: number}>} = {};
      for (const mr of info.missingResidues!) {
        if (!byChain[mr.chain]) byChain[mr.chain] = [];
        byChain[mr.chain].push({resName: mr.resName, resSeq: mr.resSeq});
      }
      const innerAcc = DG.Accordion.create();
      for (const [chain, residues] of Object.entries(byChain)) {
        innerAcc.addCountPane(`Chain ${chain}`, () => {
          residues.sort((a, b) => a.resSeq - b.resSeq);
          // Group consecutive residues into segments
          const segments: Array<Array<{resName: string; resSeq: number}>> = [];
          let current: Array<{resName: string; resSeq: number}> = [residues[0]];
          for (let i = 1; i < residues.length; i++) {
            if (residues[i].resSeq === residues[i - 1].resSeq + 1)
              current.push(residues[i]);
            else {
              segments.push(current);
              current = [residues[i]];
            }
          }
          segments.push(current);
          const parts: HTMLElement[] = [];
          for (let i = 0; i < segments.length; i++) {
            if (i > 0)
              parts.push(ui.div([], {style: {height: '8px'}}));
            parts.push(ui.divText(
              segments[i].map((r) => `${r.resName} ${r.resSeq}`).join(', '),
            ));
          }
          return ui.divV(parts);
        }, () => residues.length, false);
      }
      return innerAcc.root;
    }, () => info.missingResidues!.length, false);
  }

  const parts: HTMLElement[] = [basicTable];
  if (acc.panes.length > 0)
    parts.push(acc.root);

  return new DG.Widget(ui.divV(parts));
}
