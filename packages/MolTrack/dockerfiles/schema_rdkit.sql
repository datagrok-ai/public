-- Maintaining RDKit cartrige data in the rdk schema
create extension if not exists rdkit;
drop schema if exists rdk cascade;
create schema rdk;
drop table if exists rdk.mols;

-- rdk.mols contains ids and molecules in RDKit Mol format
select * into rdk.mols 
from (select id, mol_from_smiles(canonical_smiles::cstring) m from moltrack.compounds) tmp;

create index molidx on rdk.mols using gist(m);

-- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
-- Note that we have the "^" symbol in the body below instead of the semicolon because we
-- split the file by semicolon when executing commands (and convert back after splitting)
create or replace function insert_into_rdk_mols()
returns trigger as $$
begin
  insert into rdk.mols(id, m)
  values (
    new.id,
    public.mol_from_smiles(new.canonical_smiles::cstring)
  )^
  return null^
end^
$$ language plpgsql;

-- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
create trigger compounds_after_insert
after insert on moltrack.compounds
for each row
execute function insert_into_rdk_mols();

GRANT ALL PRIVILEGES ON SCHEMA rdk TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA rdk TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA rdk TO CURRENT_USER;