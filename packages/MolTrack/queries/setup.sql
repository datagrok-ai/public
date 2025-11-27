--name: insertAdminUser
--connection: moltrack
insert into moltrack.users (
  id, email, first_name, last_name,
  has_password, is_active, is_service_account,
  created_by, updated_by
) values (
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  'admin' || chr(64) || 'datagrok.ai', 'Admin', 'Admin',
  true, true, true,
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d'
) on conflict (email) do nothing;
--end

--name: insertSemanticTypeSynonym
--connection: moltrack
insert into moltrack.semantic_types (name, description) 
values ('Synonym', 'A semantic type representing a synonym or alternative identifier')
on conflict (name) do nothing;
--end

--name: insertProperties
--connection: moltrack
with admin as (
  select id from moltrack.users where email = 'admin' || chr(64) || 'datagrok.ai'
),
stype as (
  select id from moltrack.semantic_types where name = 'Synonym'
)
insert into moltrack.properties (created_by, updated_by, name, description, value_type, semantic_type_id, property_class, entity_type, pattern, friendly_name)
values (
  (select id from admin),
  (select id from admin),
  'corporate_compound_id', 'Official institution synonym for compounds',
  'string',
  (select id from stype),
  'DECLARED', 'COMPOUND', 'DG-{:06d}',
  'Grok ID'
), (
  (select id from admin),
  (select id from admin),
  'corporate_batch_id', 'Official institution synonym for batches',
  'string',
  (select id from stype),
  'DECLARED', 'BATCH',
  'DGB-{:06d}',
  'Grok Batch ID'
)
on conflict(name, entity_type) do nothing;
--end

--name: insertSettings
--connection: moltrack
insert into moltrack.settings (name, value, description)
values ('Compound Matching Rule',
        'ALL_LAYERS',
        'Defines the rule for matching compounds. Possible values: ALL_LAYERS (default), STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS');
--end

--name: insertStandartization
--connection: moltrack
insert into moltrack.settings (name, value, description)
values (
  'Molecule standardization rules',
  'operations:
  - type: "Cleanup"
    description: "Basic cleanup of the molecule, including removal of explicit hydrogens"
    enable: true

  - type: "FragmentParent"
    description: "Retains the largest parent fragment of the molecule."
    enable: true

  - type: "Uncharger"
    description: "Neutralizes charges on the molecule."
    enable: true
  ',
  'Defines the molecule standardization pipeline'
);
--end

--name: checkDBInitialized
--connection: moltrack
select 
	case 
		when exists (select 1 from moltrack.users where email = 'admin' || chr(64) || 'datagrok.ai')
			and exists (select 1 from moltrack.semantic_types where name = 'Synonym')
			and exists (select 1 from moltrack.properties where name = 'corporate_compound_id' and entity_type = 'COMPOUND')
			and exists (select 1 from moltrack.properties where name = 'corporate_batch_id' and entity_type = 'BATCH')
			and exists (select 1 from moltrack.settings where name = 'Compound Matching Rule')
			and exists (select 1 from moltrack.settings where name = 'Molecule standardization rules')
		then true
		else false
	end as db_initialized;
--end
