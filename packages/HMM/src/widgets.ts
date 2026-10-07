/** Context-panel widgets for a Macromolecule cell: ANARCI chain/species/
 * germlines/CDRs and Pfam domains (built-in library, gathering cutoffs). */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HmmerPool} from './pool';
import {cellResidues} from './sequences';
import {CHAIN_NAMES} from './numbering';
import {loadBuiltIn} from './models';

/** IMGT CDR definitions (all chain types). */
const IMGT_CDRS: [string, number, number][] = [['CDR1', 27, 38], ['CDR2', 56, 65], ['CDR3', 105, 117]];

function asyncWidget(fill: (root: HTMLElement) => Promise<void>): DG.Widget {
  const root = ui.div([ui.loader()]);
  fill(root).catch((e) => {
    root.replaceChildren(ui.divText(`${e instanceof Error ? e.message : e}`, {style: {color: 'var(--red-3)'}}));
  });
  return new DG.Widget(root);
}

export function anarciWidget(pool: HmmerPool, value: DG.SemanticValue): DG.Widget {
  return asyncWidget(async (root) => {
    const column = value.cell.column as DG.Column<string>;
    const {residues} = await cellResidues(column, value.cell.rowIndex);
    const [result] = await pool.anarci([['0', residues]], {scheme: 'imgt', assignGermline: true}, 1);
    if (!result.numbered || !result.details) {
      root.replaceChildren(ui.divText(result.error ? `Not numbered: ${result.error}` :
        'No antibody or T-cell receptor domain found'));
      return;
    }
    const sections: HTMLElement[] = [];
    result.numbered.forEach(([numbering], d) => {
      const details = result.details![d];
      const germlines = details.germlines && 'v_gene' in details.germlines ? details.germlines : null;
      const gene = (g: [[string, string], number] | [null, null] | undefined) =>
        g && g[0] ? `${g[0][1]} (${(g[1] * 100).toFixed(0)}%)` : '';
      const map: Record<string, string> = {
        'Chain': CHAIN_NAMES[details.chain_type] ?? details.chain_type,
        'Species': details.species,
        'V gene': gene(germlines?.v_gene),
        'J gene': gene(germlines?.j_gene),
        'E-value': details.evalue.toExponential(1),
        'Bit score': details.bitscore.toFixed(1),
      };
      for (const [name, from, to] of IMGT_CDRS) {
        map[`${name} (IMGT)`] = numbering.filter(([[n], aa]) => n >= from && n <= to && aa !== '-')
          .map(([, aa]) => aa).join('');
      }
      if (result.numbered!.length > 1) sections.push(ui.h3(`Domain ${d + 1}`));
      sections.push(ui.tableFromMap(map));
    });
    root.replaceChildren(...sections);
  });
}

export function domainsWidget(pool: HmmerPool, pkg: DG.Package, value: DG.SemanticValue): DG.Widget {
  return asyncWidget(async (root) => {
    const column = value.cell.column as DG.Column<string>;
    const {residues} = await cellResidues(column, value.cell.rowIndex);
    const library = await loadBuiltIn(pool, pkg);
    const [result] = await pool.scan(library.key, [{name: '0', residues}], {cutoff: 'ga'}, 1);
    const rows: (string | number)[][] = [];
    for (const hit of result.hits) {
      const m = library.info.models[hit.target];
      for (const d of hit.domains.filter((x) => x.reported))
        rows.push([m.name, m.accession, `${d.envFrom}-${d.envTo}`, d.bitscore.toFixed(1), d.iEvalue.toExponential(1)]);
    }
    rows.sort((a, b) => parseInt(a[2] as string) - parseInt(b[2] as string));
    root.replaceChildren(rows.length === 0 ? ui.divText('No Pfam domains (built-in library, gathering cutoffs)') :
      ui.table(rows, (r) => r, ['Domain', 'Pfam', 'Region', 'Bits', 'E-value']));
  });
}
