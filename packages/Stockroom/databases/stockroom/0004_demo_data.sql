-- A small realistic stockroom: two sites, one vendor, twenty common lab reagents with their GHS
-- hazard statements, a few containers on the shelf, and one draft purchase order. Ids are
-- md5-derived from natural keys; the auto-number counters are moved past the seeded numbers the
-- way the engine's own seed does (see DomainDdlGenerator.autoNumberSeed).

INSERT INTO stockroom.location (id, name, kind, parent_id, site)
SELECT md5('stockroom.location:' || v.site || '/' || v.name)::uuid, v.name, v.kind,
       CASE WHEN v.parent IS NULL THEN NULL ELSE md5('stockroom.location:' || v.site || '/' || v.parent)::uuid END, v.site
FROM (VALUES
  ('Main campus', 'Main campus', 'site', NULL),
  ('Main campus', 'Building A', 'building', 'Main campus'),
  ('Main campus', 'Lab 101', 'room', 'Building A'),
  ('Main campus', 'Flammables cabinet', 'cabinet', 'Lab 101'),
  ('Main campus', 'Acids cabinet', 'cabinet', 'Lab 101'),
  ('Main campus', 'Cold room', 'room', 'Building A'),
  ('Pilot plant', 'Pilot plant', 'site', NULL),
  ('Pilot plant', 'Warehouse', 'room', 'Pilot plant')
) AS v (site, name, kind, parent)
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.vendor (id, name, website, contact_email, phone) VALUES
  (md5('stockroom.vendor:Sigma-Aldrich')::uuid, 'Sigma-Aldrich', 'https://www.sigmaaldrich.com', NULL, '+1 800 325 3010')
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.substance (id, name, cas, smiles, molecular_formula)
SELECT md5('stockroom.substance:' || v.cas)::uuid, v.name, v.cas, v.smiles, v.formula
FROM (VALUES
  ('Acetone', '67-64-1', 'CC(C)=O', 'C3H6O'),
  ('Ethanol', '64-17-5', 'CCO', 'C2H6O'),
  ('Methanol', '67-56-1', 'CO', 'CH4O'),
  ('2-Propanol', '67-63-0', 'CC(C)O', 'C3H8O'),
  ('Toluene', '108-88-3', 'Cc1ccccc1', 'C7H8'),
  ('n-Hexane', '110-54-3', 'CCCCCC', 'C6H14'),
  ('Dichloromethane', '75-09-2', 'ClCCl', 'CH2Cl2'),
  ('Chloroform', '67-66-3', 'ClC(Cl)Cl', 'CHCl3'),
  ('Tetrahydrofuran', '109-99-9', 'C1CCOC1', 'C4H8O'),
  ('Ethyl acetate', '141-78-6', 'CCOC(C)=O', 'C4H8O2'),
  ('Acetonitrile', '75-05-8', 'CC#N', 'C2H3N'),
  ('Dimethyl sulfoxide', '67-68-5', 'CS(C)=O', 'C2H6OS'),
  ('Acetic acid', '64-19-7', 'CC(=O)O', 'C2H4O2'),
  ('Hydrochloric acid', '7647-01-0', 'Cl', 'HCl'),
  ('Sulfuric acid', '7664-93-9', 'OS(=O)(=O)O', 'H2SO4'),
  ('Sodium hydroxide', '1310-73-2', '[Na+].[OH-]', 'NaOH'),
  ('Sodium chloride', '7647-14-5', '[Na+].[Cl-]', 'NaCl'),
  ('Hydrogen peroxide', '7722-84-1', 'OO', 'H2O2'),
  ('Sodium azide', '26628-22-8', '[N-]=[N+]=[N-].[Na+]', 'NaN3'),
  ('Ammonium hydroxide', '1336-21-6', '[NH4+].[OH-]', 'NH4OH')
) AS v (name, cas, smiles, formula)
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.substance_hazard (id, substance_id, h_statement_id)
SELECT md5('stockroom.substance_hazard:' || v.cas || '/' || v.code)::uuid,
       md5('stockroom.substance:' || v.cas)::uuid, md5('stockroom.h_statement:' || v.code)::uuid
FROM (VALUES
  ('67-64-1', 'H225'), ('67-64-1', 'H319'), ('67-64-1', 'H336'),
  ('64-17-5', 'H225'), ('64-17-5', 'H319'),
  ('67-56-1', 'H225'), ('67-56-1', 'H301'), ('67-56-1', 'H311'), ('67-56-1', 'H331'), ('67-56-1', 'H370'),
  ('67-63-0', 'H225'), ('67-63-0', 'H319'), ('67-63-0', 'H336'),
  ('108-88-3', 'H225'), ('108-88-3', 'H304'), ('108-88-3', 'H315'), ('108-88-3', 'H336'), ('108-88-3', 'H361'), ('108-88-3', 'H373'),
  ('110-54-3', 'H225'), ('110-54-3', 'H304'), ('110-54-3', 'H315'), ('110-54-3', 'H336'), ('110-54-3', 'H361'), ('110-54-3', 'H373'), ('110-54-3', 'H411'),
  ('75-09-2', 'H315'), ('75-09-2', 'H319'), ('75-09-2', 'H336'), ('75-09-2', 'H351'),
  ('67-66-3', 'H302'), ('67-66-3', 'H315'), ('67-66-3', 'H319'), ('67-66-3', 'H331'), ('67-66-3', 'H336'), ('67-66-3', 'H351'), ('67-66-3', 'H361'), ('67-66-3', 'H372'),
  ('109-99-9', 'H225'), ('109-99-9', 'H302'), ('109-99-9', 'H319'), ('109-99-9', 'H335'), ('109-99-9', 'H351'),
  ('141-78-6', 'H225'), ('141-78-6', 'H319'), ('141-78-6', 'H336'),
  ('75-05-8', 'H225'), ('75-05-8', 'H302'), ('75-05-8', 'H312'), ('75-05-8', 'H319'), ('75-05-8', 'H332'),
  ('64-19-7', 'H226'), ('64-19-7', 'H314'),
  ('7647-01-0', 'H290'), ('7647-01-0', 'H314'), ('7647-01-0', 'H335'),
  ('7664-93-9', 'H290'), ('7664-93-9', 'H314'),
  ('1310-73-2', 'H290'), ('1310-73-2', 'H314'),
  ('7722-84-1', 'H302'), ('7722-84-1', 'H318'),
  ('26628-22-8', 'H300'), ('26628-22-8', 'H310'), ('26628-22-8', 'H373'), ('26628-22-8', 'H400'), ('26628-22-8', 'H410'),
  ('1336-21-6', 'H290'), ('1336-21-6', 'H314'), ('1336-21-6', 'H335'), ('1336-21-6', 'H400')
) AS v (cas, code)
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.container (id, label, substance_id, lot, quantity, initial_quantity, unit, received, expires, site, location_id)
SELECT md5('stockroom.container:' || v.label)::uuid, v.label, md5('stockroom.substance:' || v.cas)::uuid, v.lot,
       v.quantity, v.initial, v.unit, v.received::timestamp, v.expires::timestamp,
       'Main campus', md5('stockroom.location:Main campus/' || v.location)::uuid
FROM (VALUES
  (1000, '67-64-1', 'STBJ4821', 1.8, 2.5, 'L', '2026-01-12', '2028-01-12', 'Flammables cabinet'),
  (1001, '64-17-5', 'SHBP7730', 1.0, 1.0, 'L', '2026-02-03', '2029-02-03', 'Flammables cabinet'),
  (1002, '67-56-1', 'STBK1156', 3.2, 4.0, 'L', '2025-11-20', '2027-11-20', 'Flammables cabinet'),
  (1003, '108-88-3', 'SHBQ0304', 0.4, 1.0, 'L', '2025-09-08', '2027-09-08', 'Flammables cabinet'),
  (1004, '7664-93-9', 'STBH9962', 0.9, 1.0, 'L', '2026-03-15', '2031-03-15', 'Acids cabinet'),
  (1005, '1310-73-2', 'MKCN2277', 350, 500, 'g', '2025-06-30', '2030-06-30', 'Acids cabinet'),
  (1006, '7722-84-1', 'SZBF3390', 0.5, 1.0, 'L', '2026-04-02', '2026-10-02', 'Cold room')
) AS v (label, cas, lot, quantity, initial, unit, received, expires, location)
ON CONFLICT (id) DO NOTHING;

INSERT INTO public.domain_counters (table_id, scope_id, last)
SELECT (SELECT t.id FROM domain_tables t JOIN domain_schemas s ON s.id = t.schema_id WHERE s.name = 'stockroom' AND t.name = 'container'),
       '00000000-0000-0000-0000-000000000000'::uuid, max(label) FROM stockroom.container WHERE label IS NOT NULL
ON CONFLICT (table_id, scope_id) DO UPDATE SET last = GREATEST(domain_counters.last, excluded.last);

INSERT INTO stockroom.purchase_order (id, number, vendor_id, status, ordered_on, notes) VALUES
  (md5('stockroom.purchase_order:1')::uuid, 1, md5('stockroom.vendor:Sigma-Aldrich')::uuid, 'draft', '2026-09-01'::timestamp, 'Solvent restock for Lab 101')
ON CONFLICT (id) DO NOTHING;

INSERT INTO stockroom.order_line (id, purchase_order_id, substance_id, quantity, unit, unit_price)
SELECT md5('stockroom.order_line:1/' || v.cas)::uuid, md5('stockroom.purchase_order:1')::uuid,
       md5('stockroom.substance:' || v.cas)::uuid, v.quantity, v.unit, v.price
FROM (VALUES
  ('67-64-1', 2.5, 'L', 48.90),
  ('108-88-3', 1.0, 'L', 61.20),
  ('141-78-6', 2.5, 'L', 55.40)
) AS v (cas, quantity, unit, price)
ON CONFLICT (id) DO NOTHING;

INSERT INTO public.domain_counters (table_id, scope_id, last)
SELECT (SELECT t.id FROM domain_tables t JOIN domain_schemas s ON s.id = t.schema_id WHERE s.name = 'stockroom' AND t.name = 'purchase_order'),
       '00000000-0000-0000-0000-000000000000'::uuid, max(number) FROM stockroom.purchase_order WHERE number IS NOT NULL
ON CONFLICT (table_id, scope_id) DO UPDATE SET last = GREATEST(domain_counters.last, excluded.last);
