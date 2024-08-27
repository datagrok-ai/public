import {DbColumn, DbSchema, DbTable} from "../cruddy/table";
import {DbEntityType} from "../cruddy/entity";
import {CruddyApp, CruddyConfig, CruddyEntityView} from "../cruddy/app";

const compoundStructuresTable = new DbTable({
  name: 'compound_structures',
  columns: [
    new DbColumn({name: 'standard_inchi_key', type: 'string'}),
    new DbColumn({name: 'molregno', type: 'bigint'}),
    new DbColumn({name: 'molfile', type: 'string'}),
    new DbColumn({name: 'standard_inchi', type: 'string'}),
    new DbColumn({name: 'canonical_smiles', type: 'string'})
  ]
});

const researchCompaniesTable = new DbTable({
  name: 'research_companies',
  columns: [
    new DbColumn({name: 'company', type: 'string'}),
    new DbColumn({name: 'country', type: 'string'}),
    new DbColumn({name: 'previous_company', type: 'string'}),
    new DbColumn({name: 'co_stem_id', type: 'int'}),
    new DbColumn({name: 'res_stem_id', type: 'int'})
  ]
});

const drugMechanismTable = new DbTable({
  name: 'drug_mechanism',
  columns: [
    new DbColumn({name: 'mechanism_of_action', type: 'string'}),
    new DbColumn({name: 'disease_efficacy', type: 'int'}),
    new DbColumn({name: 'mechanism_comment', type: 'string'}),
    new DbColumn({name: 'selectivity_comment', type: 'string'}),
    new DbColumn({name: 'binding_site_comment', type: 'string'}),
    new DbColumn({name: 'variant_id', type: 'int'}),
    new DbColumn({name: 'tid', type: 'int'}),
    new DbColumn({name: 'site_id', type: 'int'}),
    new DbColumn({name: 'action_type', type: 'string'}),
    new DbColumn({name: 'direct_interaction', type: 'int'}),
    new DbColumn({name: 'molecular_mechanism', type: 'int'}),
    new DbColumn({name: 'mec_id', type: 'int'}),
    new DbColumn({name: 'record_id', type: 'int'}),
    new DbColumn({name: 'molregno', type: 'bigint', ref: 'compound_structures.molregno'})
  ]
});

const moleculeSynonymsTable = new DbTable({
  name: 'molecule_synonyms',
  columns: [
    new DbColumn({name: 'molregno', type: 'int', ref: 'compound_structures.molregno'}),
    new DbColumn({name: 'syn_type', type: 'string'}),
    new DbColumn({name: 'molsyn_id', type: 'int'}),
    new DbColumn({name: 'res_stem_id', type: 'int', ref: 'research_companies.res_stem_id'}),
    new DbColumn({name: 'synonyms', type: 'string'})
  ]
});

const moleculeDictionaryTable = new DbTable({
  name: 'molecule_dictionary',
  columns: [
    new DbColumn({name: 'chemical_probe', type: 'int'}),
    new DbColumn({name: 'molregno', type: 'int', ref: 'compound_structures.molregno'}),
    new DbColumn({name: 'pref_name', type: 'string'}),
    new DbColumn({name: 'chembl_id', type: 'string'}),
    new DbColumn({name: 'max_phase', type: 'double'}),
    new DbColumn({name: 'therapeutic_flag', type: 'int'}),
    new DbColumn({name: 'dosed_ingredient', type: 'int'}),
    new DbColumn({name: 'structure_type', type: 'string'}),
    new DbColumn({name: 'chebi_par_id', type: 'int'}),
    new DbColumn({name: 'molecule_type', type: 'string'}),
    new DbColumn({name: 'first_approval', type: 'int'}),
    new DbColumn({name: 'oral', type: 'int'}),
    new DbColumn({name: 'parenteral', type: 'int'}),
    new DbColumn({name: 'topical', type: 'int'}),
    new DbColumn({name: 'black_box_warning', type: 'int'}),
    new DbColumn({name: 'natural_product', type: 'int'}),
    new DbColumn({name: 'first_in_class', type: 'int'}),
    new DbColumn({name: 'chirality', type: 'int'}),
    new DbColumn({name: 'prodrug', type: 'int'}),
    new DbColumn({name: 'inorganic_flag', type: 'int'}),
    new DbColumn({name: 'usan_year', type: 'int'}),
    new DbColumn({name: 'availability_type', type: 'int'}),
    new DbColumn({name: 'usan_stem', type: 'string'}),
    new DbColumn({name: 'polymer_flag', type: 'int'}),
    new DbColumn({name: 'usan_substem', type: 'string'}),
    new DbColumn({name: 'usan_stem_definition', type: 'string'}),
    new DbColumn({name: 'indication_class', type: 'string'}),
    new DbColumn({name: 'withdrawn_flag', type: 'int'}),
    new DbColumn({name: 'orphan', type: 'int'})
  ]
});

const chemblSchema: DbSchema = new DbSchema('chembl',
  [compoundStructuresTable, researchCompaniesTable, drugMechanismTable, moleculeSynonymsTable, moleculeDictionaryTable]);

export const chemblConfig = new CruddyConfig({
  connection: 'Chembl:Chembl',
  schema: chemblSchema,
  entityTypes: [
    new DbEntityType({ type: 'Compound', table: compoundStructuresTable,
      gridColumnsNames: [
        'canonical_smiles', 'molecule_dictionary.chembl_id', 'drug_mechanism.action_type', 'drug_mechanism.mechanism_of_action'
      ],
      defaultView: 'cards',
      //searchColumns: ['drug_mechanism.action_type'],
      filters: [
        { type: 'distinct', column: 'molecule_dictionary.molecule_type' },
        { type: 'distinct', column: 'drug_mechanism.action_type' },
    ]}),
  ]
});


export const northwindApp = new CruddyApp(
  chemblConfig, [
    //new CruddyEntityView()
  ]
);