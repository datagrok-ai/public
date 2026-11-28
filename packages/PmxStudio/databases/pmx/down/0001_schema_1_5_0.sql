-- Domain schema pmx: 1.5.0 -> 1.4.1 (down)
-- ems-migration: from=97dfbdf3a077c87cf9a0e9fcdc36853a7e2243b01b44ea024eb6e7ddbed1d3c8 to=054b37bde4f31469d4e9151dd7a3245ac9e0a40189939c67c9654c6ce40786bc
-- Written by hand (see ../0001_schema_1_5_0.sql). Not lossless: owner goes back as the user's login (1.4.1 kept the
-- document's free-text owner there, which 0.1.2 reads only to list and search), copy_owner as the original's login,
-- run_by and copy_file_owner are dropped (run_by stays in the documents as runBy, which 0.1.2 ignores), the project
-- column comes back empty.

ALTER TABLE pmx.plan ALTER COLUMN document TYPE text USING document::text;
ALTER TABLE pmx.study ALTER COLUMN document TYPE text USING document::text;

ALTER TABLE pmx.plan ALTER COLUMN owner TYPE text USING owner::text;
ALTER TABLE pmx.study ALTER COLUMN owner TYPE text USING owner::text;
UPDATE pmx.plan p SET owner = u.login FROM public.users u WHERE u.id::text = p.owner;
UPDATE pmx.study s SET owner = u.login FROM public.users u WHERE u.id::text = s.owner;

ALTER TABLE pmx.plan ALTER COLUMN copy_owner TYPE text USING copy_owner::text;
ALTER TABLE pmx.study ALTER COLUMN copy_owner TYPE text USING copy_owner::text;
UPDATE pmx.plan p SET copy_owner = u.login FROM public.users u WHERE u.id::text = p.copy_owner;
UPDATE pmx.study s SET copy_owner = u.login FROM public.users u WHERE u.id::text = s.copy_owner;
UPDATE pmx.plan SET copy_owner = copy_file_owner WHERE copy_owner IS NULL AND copy_file_owner IS NOT NULL;
UPDATE pmx.study SET copy_owner = copy_file_owner WHERE copy_owner IS NULL AND copy_file_owner IS NOT NULL;
ALTER TABLE pmx.plan DROP COLUMN copy_file_owner;
ALTER TABLE pmx.study DROP COLUMN copy_file_owner;
ALTER TABLE pmx.plan DROP COLUMN run_by;
ALTER TABLE pmx.study DROP COLUMN run_by;

ALTER TABLE pmx.plan ALTER COLUMN start TYPE text USING to_char(start AT TIME ZONE 'UTC', 'YYYY-MM-DD');
ALTER TABLE pmx.study ALTER COLUMN start TYPE text USING to_char(start AT TIME ZONE 'UTC', 'YYYY-MM-DD');

ALTER TABLE pmx.plan ALTER COLUMN end_day TYPE double precision;
ALTER TABLE pmx.study ALTER COLUMN end_day TYPE double precision;
ALTER TABLE pmx.study ALTER COLUMN last_day TYPE double precision;
ALTER TABLE pmx.study ALTER COLUMN next_day TYPE double precision;

ALTER TABLE pmx.plan ADD COLUMN project text;
ALTER TABLE pmx.study ADD COLUMN project text;
