--name: chemblIdToSmiles
--friendlyName: Converters | ChEMBL to SMILES
--meta.role: converter
--meta.inputRegexp: (CHEMBL[0-9]+)
--connection: ChemblSql
--input: string id = "CHEMBL1185" { semType: CHEMBL_ID }
--output: string smiles { semType: Molecule }
--tags: unit-test
select canonical_smiles from compound_structures s
join molecule_dictionary d on s.molregno = d.molregno
where d.chembl_id = @id
--end


--name: molregnoToSmiles
--friendlyName: Converters | Molregno to SMILES
--connection: ChemblSql
--input: int molregno
--output: string smiles { semType: Molecule }
select canonical_smiles from compound_structures where molregno = @molregno
--end


--name: nameToSmiles
--friendlyName: Converters | Name to SMILES
--meta.role: converter
--meta.inputRegexp: (^[A-Z, a-z]+|\s+$)
--connection: ChemblSql
--input: string compoundName
--output: string smiles { semType: Molecule }
select cs.canonical_smiles
from compound_structures cs
inner join
(
  select compound_name, min(molregno) as molregno
  from compound_records
  group by compound_name
  having compound_name = @compoundName
) as cr on cr.molregno = cs.molregno
--end


--name: namesToSmiles
--friendlyName: Converters | Names to SMILES
--connection: ChemblSql
--input: list<string> names
WITH names AS (
    SELECT unnest as name FROM unnest(@names)
)
select canonical_smiles from names t1 left join
(
  select name, min(canonical_smiles) as canonical_smiles 
  from public.compound_records cr
  right join names on name = compound_name
  left join public.compound_structures cs on cr.molregno = cs.molregno
  group by name
) t2
on t1.name = t2.name
--end


--name: inchiKeyToChembl
--friendlyName: Converters | Inchi Key to ChEMBL
--connection: ChemblSql
--input: dataframe ids
select
  molecule_dictionary.chembl_id
from
  molecule_dictionary
  left join compound_structures on molecule_dictionary.molregno = compound_structures.molregno
  right join ids on compound_structures.standard_inchi_key = ids.key
  order by ids.order
--end


--name: inchiKeyToSmiles
--friendlyName: Converters | Inchi Key to SMILES
--connection: ChemblSql
--input: dataframe ids
select
  compound_structures.canonical_smiles
from
  compound_structures
  right join ids on compound_structures.standard_inchi_key = ids.key
  order by ids.order
--end


--name: inchiKeyToInchi
--friendlyName: Converters | Inchi Key to Inchi
--connection: ChemblSql
--input: dataframe ids
select
  compound_structures.standard_inchi
from
  compound_structures
  right join ids on compound_structures.standard_inchi_key = ids.key
  order by ids.order
--end

--name: chemblToSmiles
--friendlyName: Converters | ChEMBL to SMILES
--connection: ChemblSql
--input: dataframe ids
select
  compound_structures.canonical_smiles
from
  compound_structures
  left join molecule_dictionary on molecule_dictionary.molregno = compound_structures.molregno
  right join ids on molecule_dictionary.chembl_id = ids.key
  order by ids.order
--end


--name: chemblToInchi
--friendlyName: Converters | ChEMBL to Inchi
--connection: ChemblSql
--input: dataframe ids
select
  compound_structures.standard_inchi
from
  compound_structures
  left join molecule_dictionary on molecule_dictionary.molregno = compound_structures.molregno
  right join ids on molecule_dictionary.chembl_id = ids.key
  order by ids.order
--end


--name: chemblToInchiKey
--friendlyName: Converters | ChEMBL to Inchi Key
--connection: ChemblSql
--input: dataframe ids
select
  compound_structures.standard_inchi_key
from
  compound_structures
  left join molecule_dictionary on molecule_dictionary.molregno = compound_structures.molregno
  right join ids on molecule_dictionary.chembl_id = ids.key
  order by ids.order
--end
