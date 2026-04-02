import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RcsbGraphQLAdapter} from './rcsb-gql-adapter';
import {parsePdbHeaders, PdbHeaderInfo} from './pdb-helper';

/** Fetch PDB file from RCSB and parse its headers. */
async function fetchPdbFileHeaders(pdbId: string): Promise<PdbHeaderInfo | undefined> {
  try {
    const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) return undefined;
    const pdbText = await resp.text();
    return parsePdbHeaders(pdbText);
  } catch {
    return undefined;
  }
}

/** Fetch SMILES for compound IDs from RCSB Chemical Component Dictionary. */
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

/** Fetch abstract text from PubMed using the E-utilities API. */
async function fetchPubMedAbstract(pmid: number): Promise<string | undefined> {
  try {
    const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${pmid}&rettype=abstract&retmode=xml`;
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) return undefined;
    const xml = await resp.text();
    const match = xml.match(/<AbstractText[^>]*>([\s\S]*?)<\/AbstractText>/);
    return match ? match[1].replace(/<[^>]+>/g, '').trim() : undefined;
  } catch {
    return undefined;
  }
}


export async function pdbInfoWidget(pdbId: string): Promise<DG.Widget> {
  try {
    // Fetch GraphQL data and PDB file in parallel
    const [info, pdbInfo] = await Promise.all([
      RcsbGraphQLAdapter.getEntryInfo(pdbId),
      fetchPdbFileHeaders(pdbId).catch(() => undefined),
    ]);

    // -- Basic info table (same order as pdb panel) --
    const basicMap: {[key: string]: any} = {};
    if (info.title)
      basicMap['Description'] = info.title;
    if (info.classification)
      basicMap['Classification'] = info.classification;
    if (info.experimentalMethod)
      basicMap['Method'] = info.experimentalMethod;
    if (info.resolution != null)
      basicMap['Resolution'] = `${info.resolution} \u00C5`;
    basicMap['PDB URL'] = ui.link(
      `rcsb.org/structure/${pdbId}`,
      () => window.open(`https://www.rcsb.org/structure/${pdbId}`),
    );
    if (info.organisms)
      basicMap['Organism(s)'] = info.organisms.join(', ');
    if (info.expressionSystems)
      basicMap['Expression System'] = info.expressionSystems.join(', ');
    if (pdbInfo?.atomCount != null)
      basicMap['Atom Count'] = pdbInfo.atomCount.toString();
    if (pdbInfo?.chains)
      basicMap['Chains'] = pdbInfo.chains.join(', ');

    const basicTable = ui.divV([ui.h2('General'), ui.tableFromMap(basicMap)]);

    // -- Accordion for subsections --
    const acc = DG.Accordion.create();

    // -- Structure Details --
    const hasDetails = info.molecularWeight != null || info.atomCount != null ||
      info.disulfideBondCount != null || info.depositDate || info.releaseDate;
    if (hasDetails) {
      acc.addPane('Structure Details', () => {
        const map: {[key: string]: any} = {};
        if (info.molecularWeight != null)
          map['Molecular Weight'] = `${info.molecularWeight.toFixed(2)} kDa`;
        if (info.atomCount != null)
          map['Atom Count'] = info.atomCount.toString();
        if (info.disulfideBondCount != null)
          map['Disulfide Bonds'] = info.disulfideBondCount.toString();
        if (info.depositDate)
          map['Deposit Date'] = info.depositDate;
        if (info.releaseDate)
          map['Release Date'] = info.releaseDate;
        if (info.revisionDate)
          map['Revision Date'] = info.revisionDate;
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
        if (info.citation!.pubmedId) {
          citMap['PubMed'] = ui.link(
            info.citation!.pubmedId.toString(),
            () => window.open(`https://pubmed.ncbi.nlm.nih.gov/${info.citation!.pubmedId}`),
          );
        }
        const container = ui.divV([ui.tableFromMap(citMap)]);

        // Fetch abstract asynchronously if PubMed ID is available
        if (info.citation!.pubmedId) {
          const abstractHost = ui.div([ui.loader()]);
          container.appendChild(ui.divV([ui.label('Abstract'), abstractHost]));
          fetchPubMedAbstract(info.citation!.pubmedId).then((abstractText) => {
            abstractHost.innerHTML = '';
            if (abstractText)
              abstractHost.appendChild(ui.divText(abstractText));
            else
              abstractHost.appendChild(ui.divText('Abstract not available'));
          });
        }

        return container;
      }, false);
    }

    // Build binding affinity lookup by compound ID
    const affinityByComp: {[compId: string]: Array<{type: string; value: number; unit: string}>} = {};
    for (const ba of info.bindingAffinities ?? []) {
      if (!affinityByComp[ba.compId])
        affinityByComp[ba.compId] = [];
      affinityByComp[ba.compId].push({type: ba.type, value: ba.value, unit: ba.unit});
    }

    // -- Macromolecules --
    if (info.polymerEntities && info.polymerEntities.length > 0) {
      acc.addCountPane('Macromolecules', () => {
        const innerAcc = DG.Accordion.create();
        for (const pe of info.polymerEntities!) {
          const label = `Entity ${pe.entityId}` + (pe.description ? `: ${pe.description}` : '');
          innerAcc.addPane(label, () => {
            const map: {[key: string]: any} = {};
            if (pe.chains.length > 0)
              map['Chains'] = pe.chains.join(', ');
            if (pe.sequenceLength != null)
              map['Sequence Length'] = pe.sequenceLength.toString();
            // Residue count per chain from PDB file
            if (pdbInfo?.residueCountByChain && pe.authChains.length > 0) {
              const counts = pe.authChains
                .filter((c) => pdbInfo.residueCountByChain![c] != null)
                .map((c) => `${c}: ${pdbInfo.residueCountByChain![c]}`);
              if (counts.length > 0)
                map['Residues'] = counts.join(', ');
            }
            if (pe.type)
              map['Type'] = pe.type;
            if (pe.organism)
              map['Organism'] = pe.organism;
            if (pe.mutation)
              map['Mutation(s)'] = pe.mutation;
            else
              map['Mutation(s)'] = '0';
            // Specific mutations from PDB SEQADV
            if (pdbInfo?.mutations && pe.authChains.length > 0) {
              const entityMuts = pdbInfo.mutations.filter((m) => pe.authChains.includes(m.chain));
              if (entityMuts.length > 0)
                map['Mutation Details'] = entityMuts.map((m) => m.description).join(', ');
            }
            if (pe.uniprotIds) {
              map['UniProt'] = ui.divH(pe.uniprotIds.map((id) =>
                ui.link(id, () => window.open(`https://www.uniprot.org/uniprot/${id}`)),
              ));
            }
            return ui.tableFromMap(map);
          }, false);
        }
        return innerAcc.root;
      }, () => info.polymerEntities!.length, false);
    }

    // -- Small Molecules --
    if (info.nonpolymerEntities && info.nonpolymerEntities.length > 0) {
      acc.addCountPane('Small Molecules', () => {
        const innerAcc = DG.Accordion.create();
        for (const ne of info.nonpolymerEntities!) {
          innerAcc.addPane(ne.compId, () => {
            const map: {[key: string]: any} = {};
            if (ne.chains.length > 0)
              map['Chains'] = ne.chains.join(', ');
            if (ne.name)
              map['Name'] = ne.name;
            if (ne.formula)
              map['Formula'] = ne.formula;
            if (ne.smiles)
              map['SMILES'] = ne.smiles;
            if (ne.inchiKey)
              map['InChI Key'] = ne.inchiKey;

            // Binding affinity for this ligand
            const affinities = affinityByComp[ne.compId];
            if (affinities) {
              for (const ba of affinities)
                map[ba.type] = `${ba.value} ${ba.unit}`;
            }

            const parts: HTMLElement[] = [ui.tableFromMap(map)];

            // 2D structure
            if (ne.smiles)
              parts.push(grok.chem.drawMolecule(ne.smiles, 250, 200));

            return ui.divV(parts);
          }, false);
        }
        return innerAcc.root;
      }, () => info.nonpolymerEntities!.length, false);
    }

    // -- Binding Sites (from PDB file) --
    if (pdbInfo?.sites && pdbInfo.sites.length > 0) {
      acc.addCountPane('Binding Sites', () => {
        const innerAcc = DG.Accordion.create();
        for (const site of pdbInfo.sites!) {
          innerAcc.addPane(site.name, () => {
            return ui.divText(site.residues.join(', '));
          }, false);
        }
        return innerAcc.root;
      }, () => pdbInfo.sites!.length, false);
    }

    // -- Secondary Structure (from PDB file) --
    if (pdbInfo?.secondaryStructure) {
      const ss = pdbInfo.secondaryStructure;
      acc.addPane('Secondary Structure', () => {
        const map: {[key: string]: any} = {};
        if (ss.helices.length > 0) {
          map['Helices'] = ss.helices.length.toString();
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

    // -- Modified Residues (from PDB file) --
    if (pdbInfo?.modifiedResidues && pdbInfo.modifiedResidues.length > 0) {
      acc.addCountPane('Modified Residues', () => {
        const innerAcc = DG.Accordion.create();
        for (const mod of pdbInfo.modifiedResidues!) {
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

            // Fetch 2D structure from RCSB
            const molHost = ui.div([ui.loader()]);
            modParts.push(molHost);
            fetchLigandSmiles([mod.resName]).then((smilesMap) => {
              molHost.innerHTML = '';
              const smiles = smilesMap[mod.resName];
              if (smiles)
                molHost.appendChild(grok.chem.drawMolecule(smiles, 250, 200));
            });

            return ui.divV(modParts);
          }, false);
        }
        return innerAcc.root;
      }, () => pdbInfo.modifiedResidues!.length, false);
    }

    // -- Missing Residues (from PDB file) --
    if (pdbInfo?.missingResidues && pdbInfo.missingResidues.length > 0) {
      acc.addCountPane('Missing Residues', () => {
        const byChain: {[chain: string]: Array<{resName: string; resSeq: number}>} = {};
        for (const mr of pdbInfo.missingResidues!) {
          if (!byChain[mr.chain]) byChain[mr.chain] = [];
          byChain[mr.chain].push({resName: mr.resName, resSeq: mr.resSeq});
        }
        const innerAcc = DG.Accordion.create();
        for (const [chain, residues] of Object.entries(byChain)) {
          innerAcc.addCountPane(`Chain ${chain}`, () => {
            residues.sort((a, b) => a.resSeq - b.resSeq);
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
      }, () => pdbInfo.missingResidues!.length, false);
    }

    return new DG.Widget(ui.divV([basicTable, acc.root]));
  } catch (e: any) {
    return new DG.Widget(ui.divText(`Failed to load info for ${pdbId}`));
  }
}
