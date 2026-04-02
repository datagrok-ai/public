import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  parsePdbHeaders, fetchLigandSmiles,
  buildBindingSitesPane, buildSecondaryStructurePane,
  buildModifiedResiduesPane, buildMissingResiduesPane,
} from './pdb-helper';

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

  const basicTable = ui.divV([ui.h2('General'), ui.tableFromMap(basicMap, true)]);

  // -- Accordion for subsections --
  const acc = DG.Accordion.create('pdb-file-info');

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
      return ui.tableFromMap(map, true);
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
      return ui.tableFromMap(citMap, true);
    }, false);
  }

  // -- Macromolecules --
  if (info.entities && info.entities.length > 0) {
    acc.addCountPane('Macromolecules', () => {
      const innerAcc = DG.Accordion.create('pdb-file-macromolecules');
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
          return ui.tableFromMap(map, true);
        }, false);
      }
      return innerAcc.root;
    }, () => info.entities!.length, false);
  }

  // Determine if this is an RCSB PDB (has valid 4-char PDB ID) for SMILES fetching
  const isRcsb = !!(info.pdbId && /^\d[A-Za-z0-9]{3}$/.test(info.pdbId));

  // -- Small Molecules (Ligands) --
  if (info.ligands && info.ligands.length > 0) {
    acc.addCountPane('Small Molecules', () => {
      const innerAcc = DG.Accordion.create('pdb-file-small-molecules');
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
            tableHost.appendChild(ui.tableFromMap(map, true));
            fetchLigandSmiles([lig.id]).then((smilesMap) => {
              const smiles = smilesMap[lig.id];
              if (smiles) {
                map['SMILES'] = smiles;
                tableHost.innerHTML = '';
                tableHost.appendChild(ui.tableFromMap(map, true));
                molHost.innerHTML = '';
                molHost.appendChild(grok.chem.drawMolecule(smiles, 250, 200));
              } else
                molHost.innerHTML = '';
            });
          } else
            tableHost.appendChild(ui.tableFromMap(map, true));


          return ui.divV(parts);
        }, false);
      }
      return innerAcc.root;
    }, () => info.ligands!.length, false);
  }

  // -- Binding Sites --
  if (info.sites && info.sites.length > 0) {
    acc.addCountPane('Binding Sites', () => {
      return buildBindingSitesPane(info.sites!);
    }, () => info.sites!.length, false);
  }

  // -- Secondary Structure --
  if (info.secondaryStructure) {
    acc.addPane('Secondary Structure', () => {
      return buildSecondaryStructurePane(info.secondaryStructure!);
    }, false);
  }

  // -- Modified Residues --
  if (info.modifiedResidues && info.modifiedResidues.length > 0) {
    acc.addCountPane('Modified Residues', () => {
      return buildModifiedResiduesPane(info.modifiedResidues!, isRcsb);
    }, () => info.modifiedResidues!.length, false);
  }

  // -- Missing Residues --
  if (info.missingResidues && info.missingResidues.length > 0) {
    acc.addCountPane('Missing Residues', () => {
      return buildMissingResiduesPane(info.missingResidues!);
    }, () => info.missingResidues!.length, false);
  }

  const parts: HTMLElement[] = [basicTable];
  if (acc.panes.length > 0)
    parts.push(acc.root);

  return new DG.Widget(ui.divV(parts));
}
