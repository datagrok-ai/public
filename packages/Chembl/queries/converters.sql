--name: chemblIdToSmiles
--meta.role: converter
--meta.inputRegexp: (CHEMBL[0-9]+)
--connection: ChemblSql
--input: string id = "CHEMBL1185"
--output: string smiles { semType: Molecule }
--tags: unit-test
select canonical_smiles from compound_structures s
join molecule_dictionary d on s.molregno = d.molregno
where d.chembl_id = @id
--end

--name: molregnoToSmiles
--connection: ChemblSql
--input: int molregno
--output: string smiles { semType: Molecule }
select canonical_smiles from compound_structures where molregno = @molregno
--end

--name: nameToSmiles
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


--name: inchiKeyToChembl
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
