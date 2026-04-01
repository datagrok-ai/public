import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RcsbGraphQLAdapter} from './rcsb-gql-adapter';

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
    const info = await RcsbGraphQLAdapter.getEntryInfo(pdbId);

    // -- Basic info table (Description first) --
    const basicMap: {[key: string]: any} = {};
    if (info.title)
      basicMap['Description'] = info.title;
    if (info.experimentalMethod)
      basicMap['Method'] = info.experimentalMethod;
    if (info.resolution != null)
      basicMap['Resolution'] = `${info.resolution} \u00C5`;
    basicMap['PDB URL'] = ui.link(
      `rcsb.org/structure/${pdbId}`,
      () => window.open(`https://www.rcsb.org/structure/${pdbId}`),
    );
    if (info.classification)
      basicMap['Classification'] = info.classification;
    if (info.organisms)
      basicMap['Organism(s)'] = info.organisms.join(', ');
    if (info.expressionSystems)
      basicMap['Expression System'] = info.expressionSystems.join(', ');

    const basicTable = ui.tableFromMap(basicMap);

    // -- Accordion for subsections --
    const acc = DG.Accordion.create();

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
            if (pe.type)
              map['Type'] = pe.type;
            if (pe.organism)
              map['Organism'] = pe.organism;
            if (pe.mutation)
              map['Mutation(s)'] = pe.mutation;
            else
              map['Mutation(s)'] = '0';
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
            if (ne.inchiKey)
              map['InChI Key'] = ne.inchiKey;
            const parts: HTMLElement[] = [ui.tableFromMap(map)];
            if (ne.smiles)
              parts.push(grok.chem.drawMolecule(ne.smiles, 250, 200));
            return ui.divV(parts);
          }, false);
        }
        return innerAcc.root;
      }, () => info.nonpolymerEntities!.length, false);
    }

    return new DG.Widget(ui.divV([basicTable, acc.root]));
  } catch (e: any) {
    return new DG.Widget(ui.divText(`Failed to load info for ${pdbId}`));
  }
}
