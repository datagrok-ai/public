import {Section} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainSearch} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions, DomainRibbon} from '@datagrok-libraries/ems/src/dg';
import {MolTrackEms} from './app';
import {RegisterDialog} from './register-dialog';
import {MolTrackTables} from './tables';

/** The compound page: what a chemist reads a compound by first, the registration internals folded under it; a
 * structure drawn on a new compound that is registered already shows the compound it is. Compounds arrive through
 * registration: Import runs the bulk registration, not a plain import. */
export class CompoundsApp extends DomainApp {
  static readonly FIELDS = ['molregno', 'original_structure', 'structure', 'common_name', 'cas', 'epa_compound_id',
    'usan', 'iupac_name', 'source', 'source_compound_code', 'mol_logp', 'formula', 'inchikey'];
  static readonly DETAILS = ['hash_mol', 'match_rule', 'hash_all_layers', 'hash_stereo_insensitive',
    'hash_tautomer_insensitive', 'inchi'];
  static readonly SEARCH_HINT =
    'Searches names, CAS numbers, IDs and formulas. To find a structure, add a Structure condition under Filters';

  constructor(options: DomainAppOptions) {
    super(options);
    MolTrackEms.aroundForm(this, (form, source) => {
      const registered = MolTrackEms.compoundSection('Already registered', () => {
        source.rows.items.value;
        const row = source.currentRow.value;
        return row === null ? null : MolTrackTables.duplicateOf(row)?.id ?? null;
      });
      const details = new Section({title: 'Registration details', expanded: false,
        hint: 'Registration hashes, matching rule, InChI', persistKey: 'moltrack-ems-registration-details'});
      details.add(domains.form(source, {include: CompoundsApp.DETAILS}).root);
      form.root.before(registered.root);
      form.root.after(details.root);
    });
  }

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
