import {DbColumn, DbSchema, DbTable} from "../cruddy/table";
import {DbEntityType} from "../cruddy/entity";
import {CruddyApp, CruddyConfig, CruddyEntityView} from "../cruddy/app";

// Cruddy template
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

// Cruddy template
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

// Cruddy template
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

// Cruddy template
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


const chemblSchema: DbSchema = new DbSchema('chembl',
  [compoundStructuresTable, researchCompaniesTable, drugMechanismTable, moleculeSynonymsTable]);

export const chemblConfig = new CruddyConfig({
  connection: 'Chembl:Chembl',
  schema: chemblSchema,
  entityTypes: [
    new DbEntityType({ type: 'Compound', table: compoundStructuresTable,
      gridColumnsNames: [
        'canonical_smiles', 'drug_mechanism.action_type', 'drug_mechanism.mechanism_of_action'
      ],
      filters: [
        // { type: 'distinct', column: 'shipcountry'},
        // { type: 'combo', column: 'shipcity'},
        // { type: 'range', column: 'orderid'},
        // { type: 'expression', column: 'shipcity'},
    ]}),
  ]
});


export const northwindApp = new CruddyApp(
  chemblConfig, [
    //new CruddyEntityView()
  ]
);