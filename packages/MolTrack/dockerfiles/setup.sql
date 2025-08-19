BEGIN;

INSERT INTO moltrack.users (
  id, email, first_name, last_name,
  has_password, is_active, is_service_account,
  created_by, updated_by
) VALUES (
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  'admin@datagrok.ai', 'Admin', 'Admin',
  true, true, true,
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d'
);

INSERT INTO moltrack.semantic_types (name, description) 
VALUES ('Synonym', 'A semantic type representing a synonym or alternative identifier')
ON CONFLICT (name) DO NOTHING;

with ADMIN AS (
  SELECT id FROM moltrack.users WHERE email = 'admin@datagrok.ai'
),
STYPE AS (
  SELECT id FROM moltrack.semantic_types WHERE name = 'Synonym'
)
INSERT INTO moltrack.properties (created_by, updated_by, name, description, value_type, semantic_type_id, property_class, entity_type, pattern)
VALUES (
  (SELECT id FROM ADMIN),
  (SELECT id FROM ADMIN),
  'corporate_compound_id', 'Official institution synonym for compounds',
  'string', 
  (SELECT id FROM STYPE), 
  'DECLARED', 'COMPOUND', 'DG-{:06d}'
), (
  (SELECT id FROM ADMIN),
  (SELECT id FROM ADMIN),
  'corporate_batch_id', 'Official institution synonym for batches',
  'string', 
  (SELECT id FROM STYPE), 
  'DECLARED', 'BATCH',
  'DGB-{:06d}'
);

INSERT INTO moltrack.settings (name, value, description)
VALUES ('Compound Matching Rule',
        'ALL_LAYERS',
        'Defines the rule for matching compounds. Possible values: ALL_LAYERS (default), STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS');

COMMIT;