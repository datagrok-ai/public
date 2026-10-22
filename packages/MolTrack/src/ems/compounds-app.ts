import {DomainApp, DomainSearch} from '@datagrok-libraries/ems/src/dg';
import type {DomainRibbon} from '@datagrok-libraries/ems/src/dg';

/** The compound page: what a chemist reads a compound by first, the registration identity folded under it. */
export class CompoundsApp extends DomainApp {
  static readonly FIELDS = ['molregno', 'original_structure', 'common_name', 'cas', 'epa_compound_id', 'usan',
    'iupac_name', 'source', 'source_compound_code', 'mol_logp', 'formula', 'inchikey'];
  static readonly DETAILS = ['structure', 'hash_mol', 'match_rule', 'hash_all_layers', 'hash_stereo_insensitive',
    'hash_tautomer_insensitive', 'inchi'];
  static readonly SEARCH_HINT =
    'Searches names, CAS numbers, EPA IDs and formulas. For a Grok ID or a structure, add a condition under Filters';

  ribbon(): DomainRibbon {
    const ribbon = super.ribbon();
    for (const tool of ribbon.tools) {
      if (tool instanceof DomainSearch)
        tool.tooltipText = CompoundsApp.SEARCH_HINT;
    }
    return ribbon;
  }
}
