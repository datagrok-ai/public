--name: _compoundNames
--friendlyName: Suggestions | Compound Names
--description: Provides autocomplete suggestions for compound names matching a substring pattern.
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct compound_name from compound_records
where compound_name ilike '%' || @sub || '%'
limit 50
--end

--name: _organisms
--friendlyName: Suggestions | Organisms
--description: Provides autocomplete suggestions for organism names matching a substring pattern.
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct organism from target_dictionary
where organism ilike '%' || @sub || '%'
limit 50
--end

--name: _proteinTypes
--friendlyName: Suggestions | Protein Types
--description: Provides autocomplete suggestions for protein classification short names matching a substring pattern.
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct short_name from protein_classification
where short_name ilike '%' || @sub || '%'
limit 50
--end

--name: _targetTypes
--friendlyName: Suggestions | Target Types
--description: Provides autocomplete suggestions for target types matching a substring pattern.
--connection: Chembl
--input: string sub
select target_type from target_type
where target_type ilike '%' || @sub || '%'
limit 50
--end

--name: _assayTypes
--friendlyName: Suggestions | Assay Types
--description: Provides autocomplete suggestions for assay types matching a substring pattern.
--connection: Chembl
--input: string sub
select assay_type from assay_type
where assay_type ilike '%' || @sub || '%'
limit 50
--end

--name: _relationshipTypes
--friendlyName: Suggestions | Relationship Types
--description: Provides autocomplete suggestions for relationship types matching a substring pattern.
--connection: Chembl
--input: string sub
select relationship_type from relationship_type
where relationship_type ilike '%' || @sub || '%'
limit 50
--end