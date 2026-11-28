-- Domain schema pmx: 1.4.1 -> 1.5.0
-- ems-migration: from=054b37bde4f31469d4e9151dd7a3245ac9e0a40189939c67c9654c6ce40786bc to=97dfbdf3a077c87cf9a0e9fcdc36853a7e2243b01b44ea024eb6e7ddbed1d3c8
-- Written by hand on the pattern of Chem's mpo 0001_properties_json.sql and Grit's 0001_issue_autonumber.sql: the
-- header's snapshot hashes, and the physical types (uuid, timestamptz, integer, jsonb), must be checked against the
-- scaffold the platform's schema snapshot tool gives for 1.4.1 -> 1.5.0 before the package is published anywhere
-- with 1.4.1 rows (dev). Each statement is run once, in order, in the deploy's transaction.
-- Registry-only changes (no SQL): document loses its "editor"; descriptions; the new run_by and copy_file_owner
-- columns (created by the deploy; added here too, IF NOT EXISTS, since the data moves below need them).

ALTER TABLE pmx.plan ADD COLUMN IF NOT EXISTS run_by text;
ALTER TABLE pmx.plan ADD COLUMN IF NOT EXISTS copy_file_owner text;
ALTER TABLE pmx.study ADD COLUMN IF NOT EXISTS run_by text;
ALTER TABLE pmx.study ADD COLUMN IF NOT EXISTS copy_file_owner text;

-- change-type plan.document, study.document: type string -> json [manual]
-- Every value the app stored is the JSON.stringify of a SavedPlan / SavedStudy (store.ts columnsOf). A text that is not
-- JSON object (a row written by hand, e.g. the 11-24 size probe's 'x' * n, trashed) is kept whole inside an object:
-- the deploy's precheck refuses a json column holding anything but objects, after this script has committed.
CREATE FUNCTION pg_temp.pmx_as_jsonb(t text) RETURNS jsonb LANGUAGE plpgsql IMMUTABLE AS $f$
  DECLARE v jsonb;
  BEGIN v := t::jsonb; RETURN CASE WHEN jsonb_typeof(v) = 'object' THEN v ELSE jsonb_build_object('value', v) END;
  EXCEPTION WHEN invalid_text_representation THEN RETURN jsonb_build_object('text', t); END $f$;
ALTER TABLE pmx.plan ALTER COLUMN document TYPE jsonb USING pg_temp.pmx_as_jsonb(document);
ALTER TABLE pmx.study ALTER COLUMN document TYPE jsonb USING pg_temp.pmx_as_jsonb(document);

-- change-type plan.owner, study.owner: type string -> user [manual]
-- owner held the document's free-text owner (by default the author's login). The row's user is its author (every
-- row was inserted by whoever owns it: a save, an import or a Save as my copy). A text that is not the login of any
-- user named whoever runs the study: it moves to run_by, and into the document as runBy (Study details' Run by). A
-- login never does, the author's or another's: a copy kept the original owner's login there, and a login is a user, not
-- free text.
UPDATE pmx.plan p SET run_by = p.owner,
    document = CASE WHEN p.document ? 'study' THEN jsonb_set(p.document, '{study,runBy}', to_jsonb(p.owner))
      ELSE p.document END
  WHERE coalesce(p.owner, '') <> ''
    AND NOT EXISTS (SELECT 1 FROM public.users u WHERE lower(u.login) = lower(p.owner));
UPDATE pmx.study s SET run_by = s.owner,
    document = CASE WHEN s.document ? 'meta' THEN jsonb_set(s.document, '{meta,runBy}', to_jsonb(s.owner))
      ELSE s.document END
  WHERE coalesce(s.owner, '') <> ''
    AND NOT EXISTS (SELECT 1 FROM public.users u WHERE lower(u.login) = lower(s.owner));
ALTER TABLE pmx.plan ALTER COLUMN owner TYPE uuid USING author_id;
ALTER TABLE pmx.study ALTER COLUMN owner TYPE uuid USING author_id;

-- change-type plan.copy_owner, study.copy_owner: type string -> user [manual]
-- copy_owner held a login or name (a file brought in as the user's copy, or the original's owner field) or the
-- original's author id. A user id or a login becomes the user; any other text is the owner a file named, and moves to
-- copy_file_owner.
UPDATE pmx.plan p SET copy_owner = u.id::text FROM public.users u WHERE lower(u.login) = lower(p.copy_owner);
UPDATE pmx.study s SET copy_owner = u.id::text FROM public.users u WHERE lower(u.login) = lower(s.copy_owner);
UPDATE pmx.plan SET copy_file_owner = copy_owner, copy_owner = NULL
  WHERE coalesce(copy_owner, '') <> '' AND copy_owner !~* '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$';
UPDATE pmx.study SET copy_file_owner = copy_owner, copy_owner = NULL
  WHERE coalesce(copy_owner, '') <> '' AND copy_owner !~* '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$';
ALTER TABLE pmx.plan ALTER COLUMN copy_owner TYPE uuid USING nullif(copy_owner, '')::uuid;
ALTER TABLE pmx.study ALTER COLUMN copy_owner TYPE uuid USING nullif(copy_owner, '')::uuid;

-- change-type plan.start, study.start: type string -> datetime [manual]
-- start held day 0 as YYYY-MM-DD (or nothing): midnight UTC of that day.
ALTER TABLE pmx.plan ALTER COLUMN start TYPE timestamptz
  USING CASE WHEN start ~ '^\d{4}-\d{2}-\d{2}$' THEN (start || 'T00:00:00Z')::timestamptz END;
ALTER TABLE pmx.study ALTER COLUMN start TYPE timestamptz
  USING CASE WHEN start ~ '^\d{4}-\d{2}-\d{2}$' THEN (start || 'T00:00:00Z')::timestamptz END;

-- change-type plan.end_day, study.end_day, study.last_day, study.next_day: type float -> int [manual]
-- Whole study days: a planned end the model gave with a fraction is rounded, as the app now rounds it on write.
ALTER TABLE pmx.plan ALTER COLUMN end_day TYPE integer USING round(end_day)::integer;
ALTER TABLE pmx.study ALTER COLUMN end_day TYPE integer USING round(end_day)::integer;
ALTER TABLE pmx.study ALTER COLUMN last_day TYPE integer USING round(last_day)::integer;
ALTER TABLE pmx.study ALTER COLUMN next_day TYPE integer USING round(next_day)::integer;

-- drop-column plan.project, study.project [manual]
-- Nothing is lost: the project is also in each document (study.project / meta.project), which the app reads as one
-- more tag from 0.1.3 on (study-meta.ts parseMeta); the row's tags take it at its next save.
ALTER TABLE pmx.plan DROP COLUMN project;
ALTER TABLE pmx.study DROP COLUMN project;
