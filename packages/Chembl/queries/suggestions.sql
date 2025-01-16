--name: _compoundNames
--friendlyName: Suggestions | Compound Names
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct compound_name from compound_records
where compound_name ilike '%' || @sub || '%'
limit 50
--end

--name: _organisms
--friendlyName: Suggestions | Organisms
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct organism from target_dictionary
where organism ilike '%' || @sub || '%'
limit 50
--end

--name: _proteinTypes
--friendlyName: Suggestions | Protein Types
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct short_name from protein_classification
where short_name ilike '%' || @sub || '%'
limit 50
--end

--name: _targetTypes
--friendlyName: Suggestions | Target Types
--connection: Chembl
--input: string sub
select target_type from target_type
where target_type ilike '%' || @sub || '%'
limit 50
--end

--name: _assayTypes
--friendlyName: Suggestions | Assay Types
--connection: Chembl
--input: string sub
select assay_type from assay_type
where assay_type ilike '%' || @sub || '%'
limit 50
--end

--name: _relationshipTypes
--friendlyName: Suggestions | Relationship Types
--connection: Chembl
--input: string sub
select relationship_type from relationship_type
where relationship_type ilike '%' || @sub || '%'
limit 50
--end