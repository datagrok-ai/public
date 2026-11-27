-- One-off: the platform's earlier models of package methods become Forge models, promoted as the server promotes a row.
WITH src AS (
  SELECT uuid_generate_v4() AS id, m.id::text AS legacy_id, m.author_id, g.id AS author_group,
         COALESCE(m.created_on AT TIME ZONE 'UTC', now()) AS created_on,
         COALESCE(COALESCE(m.updated_on, m.created_on) AT TIME ZONE 'UTC', now()) AS updated_on,
         COALESCE(NULLIF(trim(m.friendly_name), ''), m.name) AS name, m.description,
         substr(m.source, strpos(m.source, ': ') + 2) AS engine_name,
         left(m.source, strpos(m.source, ': ') - 1) AS engine_namespace,
         CASE WHEN e.value_string LIKE '%.Script' THEN 'script' ELSE 'function' END AS engine_kind,
         CASE WHEN o.type IN ('string', 'bool') THEN 'classification' ELSE 'regression' END AS task,
         o.name AS target_name,
         jsonb_strip_nulls(jsonb_build_object('name', o.name, 'type', o.type, 'semType', o.sem_type)) AS target,
         f.columns,
         CASE WHEN json_typeof(m.options) = 'object' THEN m.options::jsonb ELSE '{}' END AS options,
         t.name AS dataset_name
  FROM public.models m
  JOIN public.entities en ON en.id = m.id AND NOT en.is_deleted
  JOIN public.users u ON u.id = m.author_id
  JOIN public.groups g ON g.id = u.group_id
  JOIN LATERAL (SELECT c.name, c.type, c.sem_type
                FROM public.pm_output_columns p JOIN public.table_columns c ON c.id = p.column_id
                WHERE p.model_id = m.id AND c.name <> '' AND c.type <> ''
                ORDER BY p.ctid LIMIT 1) o ON true
  JOIN LATERAL (SELECT jsonb_agg(jsonb_strip_nulls(jsonb_build_object('name', d.name, 'type', d.type,
                                                                       'semType', d.sem_type)) ORDER BY d.pos) AS columns
                FROM (SELECT DISTINCT ON (c.name) c.name, c.type, c.sem_type, p.ctid AS pos
                      FROM public.pm_input_columns p JOIN public.table_columns c ON c.id = p.column_id
                      WHERE p.model_id = m.id AND c.name <> '' AND c.type <> ''
                      ORDER BY c.name, p.ctid) d) f ON f.columns IS NOT NULL
  LEFT JOIN public.tables t ON t.id = m.trained_on_id
  LEFT JOIN LATERAL (SELECT v.value_string
                     FROM public.entity_property_values v JOIN public.entity_properties p ON p.id = v.property_id
                     WHERE v.entity_id = m.id AND p.name = 'engineType' LIMIT 1) e ON true
  WHERE strpos(m.source, ': ') > 1
    AND trim(substr(m.source, strpos(m.source, ': ') + 2)) <> ''
    AND COALESCE(NULLIF(trim(m.friendly_name), ''), m.name) <> ''
    AND NOT EXISTS (SELECT 1 FROM model x WHERE x.legacy_id = m.id::text)),
ins AS (
  INSERT INTO model (id, legacy_id, author_id, created_on, updated_on, name, description, engine_name,
                     engine_namespace, engine_kind, task, target_name, target, features, feature_count, options,
                     storage_mode, dataset_name)
  SELECT id, legacy_id, author_id, created_on, updated_on, name, description, engine_name, engine_namespace,
         engine_kind, task, target_name, target, jsonb_build_object('columns', columns),
         jsonb_array_length(columns), options, 'none', dataset_name
  FROM src),
ent AS (
  INSERT INTO public.entities (id, entity_type_id, handle, name, friendly_name, is_deleted)
  SELECT src.id, t.id, 'forge.model:' || src.id, 'forge.model:' || src.id, src.name, false
  FROM src JOIN public.entity_types t ON t.name = 'forge.model'
  ON CONFLICT (id) DO NOTHING)
INSERT INTO public.permissions (id, entity_id, user_group_id, permission_id)
SELECT uuid_generate_v5('0bb1a16e-df07-45f5-99e0-8e67000ecc22',
                        'DomainGrant:' || src.id || ':' || src.author_group || ':' || p.name),
       src.id, src.author_group, p.id
FROM src
JOIN public.entity_types pt ON pt.name = 'Project'
JOIN public.entity_types_permissions p ON p.entity_type_id = pt.id AND p.name IN ('View', 'Edit', 'Delete', 'Share')
ON CONFLICT DO NOTHING;
