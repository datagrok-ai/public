import {DomainApp, DomainSearch} from '@datagrok-libraries/ems/src/dg';
import type {DomainRibbon} from '@datagrok-libraries/ems/src/dg';
import {RegisterDialog} from './register-dialog';

/** The compound page: what a chemist reads a compound by first, the registration identity folded under it.
 * Compounds arrive through registration: Import runs the bulk registration, not a plain import. */
export class CompoundsApp extends DomainApp {
  static readonly FIELDS = ['molregno', 'original_structure', 'common_name', 'cas', 'epa_compound_id', 'usan',
    'iupac_name', 'source', 'source_compound_code', 'mol_logp', 'formula', 'inchikey'];
  static readonly DETAILS = ['structure', 'hash_mol', 'match_rule', 'hash_all_layers', 'hash_stereo_insensitive',
    'hash_tautomer_insensitive', 'inchi'];
  static readonly SEARCH_HINT =
    'Searches names, CAS numbers, IDs and formulas. To find a structure, add a Structure condition under Filters';

  async importRows(): Promise<boolean> {
    return await RegisterDialog.open() !== null;
  }

  ribbon(): DomainRibbon {
    const ribbon = super.ribbon();
    for (const tool of ribbon.tools) {
      if (tool instanceof DomainSearch)
        tool.tooltipText = CompoundsApp.SEARCH_HINT;
    }
    return ribbon;
  }
}
