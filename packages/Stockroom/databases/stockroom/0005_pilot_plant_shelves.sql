-- The Pilot plant branch of the location tree, so both roots carry a real subtree: the
-- warehouse gets two shelves and two containers stand on them. What this is for is the
-- location tree driving the container list — selecting "Pilot plant" must answer rows,
-- not an empty page (`location_id under "<id>"`, schema.json "hierarchy": true on location).
-- The labels start at 1100 so a stand where someone has already made containers by hand
-- (the business key is the label) cannot collide with the seed.

INSERT INTO stockroom.location (id, name, kind, parent_id, site)
SELECT md5('stockroom.location:' || v.site || '/' || v.name)::uuid, v.name, v.kind,
       md5('stockroom.location:' || v.site || '/' || v.parent)::uuid, v.site
FROM (VALUES
  ('Pilot plant', 'Solvent store', 'cabinet', 'Warehouse'),
  ('Pilot plant', 'Reagent store', 'cabinet', 'Warehouse')
) AS v (site, name, kind, parent)
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.container (id, label, substance_id, lot, quantity, initial_quantity, unit, received, expires, site, location_id)
SELECT md5('stockroom.container:' || v.label)::uuid, v.label, md5('stockroom.substance:' || v.cas)::uuid, v.lot,
       v.quantity, v.initial, v.unit, v.received::timestamp, v.expires::timestamp,
       'Pilot plant', md5('stockroom.location:Pilot plant/' || v.location)::uuid
FROM (VALUES
  (1100, '141-78-6', 'STBM4417', 18.0, 20.0, 'L', '2026-05-18', '2028-05-18', 'Solvent store'),
  (1101, '7647-14-5', 'MKCP8831', 22.5, 25.0, 'kg', '2026-06-02', '2031-06-02', 'Reagent store')
) AS v (label, cas, lot, quantity, initial, unit, received, expires, location)
ON CONFLICT (id) DO NOTHING;

INSERT INTO public.domain_counters (table_id, scope_id, last)
SELECT (SELECT t.id FROM domain_tables t JOIN domain_schemas s ON s.id = t.schema_id WHERE s.name = 'stockroom' AND t.name = 'container'),
       '00000000-0000-0000-0000-000000000000'::uuid, max(label) FROM stockroom.container WHERE label IS NOT NULL
ON CONFLICT (table_id, scope_id) DO UPDATE SET last = GREATEST(domain_counters.last, excluded.last);
