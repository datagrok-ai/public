-- Catalog fixture for db-explorer: three working catalogs. `research` and `archive` hold the same
-- chem.compounds table with the same ids and different names, so a query that drops the catalog
-- returns visibly wrong rows; `sales` and the two hundred `tenant_NNN` catalogs stand beside them.
--
-- catalogs-db.mjs runs it as sa with `-v owner=catalogs_demo`: every database it makes is given to
-- that login at once, and it drops or changes only databases of its own names, owned by that login or
-- still by sa. A database of one of those names that another login owns stops it before anything
-- changes, so the databases other logins keep on the same server are never touched. It can be run
-- again: it re-creates its three working catalogs and adds the two hundred others only where they are
-- missing.
:on error exit
SET NOCOUNT ON;
GO
IF SUSER_ID(N'$(owner)') IS NULL
  RAISERROR(N'The login $(owner) does not exist: the fixture gives its databases to it.', 16, 1);
GO
DECLARE @mine TABLE (name SYSNAME COLLATE DATABASE_DEFAULT PRIMARY KEY);
INSERT INTO @mine VALUES (N'research'), (N'archive'), (N'sales');
DECLARE @i INT = 1;
WHILE @i <= 200
BEGIN
  INSERT INTO @mine VALUES (N'tenant_' + RIGHT('000' + CAST(@i AS VARCHAR(3)), 3));
  SET @i += 1;
END
DECLARE @others NVARCHAR(2000) = (SELECT STRING_AGG(CAST(d.name AS NVARCHAR(200)), N', ')
  FROM sys.databases d JOIN @mine m ON m.name = d.name COLLATE DATABASE_DEFAULT
  WHERE ISNULL(SUSER_SNAME(d.owner_sid), N'') NOT IN (N'$(owner)', N'sa'));
IF @others IS NOT NULL
  RAISERROR(N'%s: of the fixture''s names, but another login owns them. The fixture leaves them alone and stops.',
    16, 1, @others);
GO
IF DB_ID('research') IS NOT NULL BEGIN ALTER DATABASE research SET SINGLE_USER WITH ROLLBACK IMMEDIATE; DROP DATABASE research; END
IF DB_ID('archive') IS NOT NULL BEGIN ALTER DATABASE archive SET SINGLE_USER WITH ROLLBACK IMMEDIATE; DROP DATABASE archive; END
IF DB_ID('sales') IS NOT NULL BEGIN ALTER DATABASE sales SET SINGLE_USER WITH ROLLBACK IMMEDIATE; DROP DATABASE sales; END
GO
CREATE DATABASE research;
ALTER AUTHORIZATION ON DATABASE::research TO [$(owner)];
CREATE DATABASE archive;
ALTER AUTHORIZATION ON DATABASE::archive TO [$(owner)];
CREATE DATABASE sales;
ALTER AUTHORIZATION ON DATABASE::sales TO [$(owner)];
GO

------------------------------------------------------------------------------------------------
USE research;
GO
CREATE SCHEMA core;
GO
CREATE SCHEMA chem;
GO
CREATE SCHEMA bio;
GO
CREATE TABLE core.projects (
  project_id INT NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  lead VARCHAR(64) NOT NULL,
  status VARCHAR(16) NOT NULL
);
CREATE TABLE chem.compounds (
  compound_id VARCHAR(16) NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  smiles VARCHAR(256) NOT NULL,
  mol_weight FLOAT NOT NULL,
  project_id INT NOT NULL CONSTRAINT fk_compounds_project REFERENCES core.projects (project_id)
);
CREATE TABLE chem.batches (
  batch_id VARCHAR(16) NOT NULL PRIMARY KEY,
  compound_id VARCHAR(16) NOT NULL CONSTRAINT fk_batches_compound REFERENCES chem.compounds (compound_id),
  purity FLOAT NOT NULL,
  synthesized_on DATE NOT NULL,
  chemist VARCHAR(64) NOT NULL
);
CREATE TABLE bio.assays (
  assay_id INT NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  target VARCHAR(32) NOT NULL,
  assay_type VARCHAR(32) NOT NULL
);
CREATE TABLE bio.assay_results (
  result_id INT NOT NULL PRIMARY KEY,
  batch_id VARCHAR(16) NOT NULL CONSTRAINT fk_results_batch REFERENCES chem.batches (batch_id),
  assay_id INT NOT NULL CONSTRAINT fk_results_assay REFERENCES bio.assays (assay_id),
  ic50_nm FLOAT NOT NULL,
  measured_on DATE NOT NULL
);
GO
INSERT INTO core.projects VALUES
  (1, 'Kinase inhibitors', 'Maria Chen', 'active'),
  (2, 'CNS penetrants', 'Omar Haddad', 'active'),
  (3, 'Anti-inflammatory', 'Ingrid Solberg', 'on hold');
INSERT INTO chem.compounds VALUES
  ('DGC-000001', 'Aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 180.16, 3),
  ('DGC-000002', 'Ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 206.28, 3),
  ('DGC-000003', 'Naproxen', 'COc1ccc2cc(ccc2c1)C(C)C(=O)O', 230.26, 3),
  ('DGC-000004', 'Paracetamol', 'CC(=O)Nc1ccc(O)cc1', 151.16, 3),
  ('DGC-000005', 'Caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 194.19, 2),
  ('DGC-000006', 'Nicotine', 'CN1CCCC1c1cccnc1', 162.23, 2),
  ('DGC-000007', 'Diazepam', 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12', 284.74, 2),
  ('DGC-000008', 'Fluoxetine', 'CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1', 309.33, 2),
  ('DGC-000009', 'Imatinib', 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1', 493.6, 1),
  ('DGC-000010', 'Gefitinib', 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1', 446.9, 1),
  ('DGC-000011', 'Erlotinib', 'COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC', 393.44, 1),
  ('DGC-000012', 'Sorafenib', 'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(c3)C(F)(F)F)cc2)ccn1', 464.83, 1);
INSERT INTO chem.batches (batch_id, compound_id, purity, synthesized_on, chemist)
SELECT 'DGB-' + RIGHT('000000' + CAST((ROW_NUMBER() OVER (ORDER BY c.compound_id, n.n)) AS VARCHAR(6)), 6),
       c.compound_id,
       95.0 + ((ABS(CHECKSUM(c.compound_id)) + n.n * 7) % 49) / 10.0,
       DATEADD(DAY, (ABS(CHECKSUM(c.compound_id)) % 200) + n.n * 30, '2024-01-15'),
       CASE n.n WHEN 1 THEN 'Lena Fischer' ELSE 'Tariq Noor' END
FROM chem.compounds c CROSS JOIN (VALUES (1), (2)) n(n);
INSERT INTO bio.assays VALUES
  (1, 'EGFR enzymatic', 'EGFR', 'biochemical'),
  (2, 'ABL1 enzymatic', 'ABL1', 'biochemical'),
  (3, 'COX-2 inhibition', 'PTGS2', 'biochemical'),
  (4, 'HepG2 viability', 'cell line', 'cellular');
INSERT INTO bio.assay_results (result_id, batch_id, assay_id, ic50_nm, measured_on)
SELECT ROW_NUMBER() OVER (ORDER BY b.batch_id, a.assay_id),
       b.batch_id,
       a.assay_id,
       ROUND(1 + (ABS(CHECKSUM(b.batch_id + CAST(a.assay_id AS VARCHAR(4)))) % 99000) / 10.0, 1),
       DATEADD(DAY, a.assay_id * 3, b.synthesized_on)
FROM chem.batches b CROSS JOIN bio.assays a
WHERE (ABS(CHECKSUM(b.batch_id)) + a.assay_id) % 4 <> 0;
GO

------------------------------------------------------------------------------------------------
USE archive;
GO
CREATE SCHEMA core;
GO
CREATE SCHEMA chem;
GO
CREATE TABLE core.projects (
  project_id INT NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  lead VARCHAR(64) NOT NULL,
  status VARCHAR(16) NOT NULL
);
CREATE TABLE chem.compounds (
  compound_id VARCHAR(16) NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  smiles VARCHAR(256) NOT NULL,
  mol_weight FLOAT NOT NULL,
  project_id INT NOT NULL CONSTRAINT fk_compounds_project REFERENCES core.projects (project_id)
);
GO
INSERT INTO core.projects VALUES (1, 'Archived 2019 screen', 'Retired Lead', 'closed');
INSERT INTO chem.compounds
SELECT compound_id, 'ARCHIVED ' + name, 'C', 16.04, 1 FROM research.chem.compounds;
GO

------------------------------------------------------------------------------------------------
USE sales;
GO
CREATE TABLE dbo.customers (
  customer_id INT NOT NULL PRIMARY KEY,
  name VARCHAR(64) NOT NULL,
  country VARCHAR(32) NOT NULL
);
CREATE TABLE dbo.orders (
  order_id INT NOT NULL PRIMARY KEY,
  customer_id INT NOT NULL CONSTRAINT fk_orders_customer REFERENCES dbo.customers (customer_id),
  ordered_on DATE NOT NULL,
  total FLOAT NOT NULL
);
GO
INSERT INTO dbo.customers VALUES (1, 'Northwind Labs', 'USA'), (2, 'Helix Bio', 'Germany'), (3, 'Kyoto Pharma', 'Japan');
INSERT INTO dbo.orders VALUES
  (1001, 1, '2025-03-02', 1250.50), (1002, 1, '2025-04-11', 310.00),
  (1003, 2, '2025-04-19', 9800.00), (1004, 3, '2025-05-07', 455.25);
GO
-- 200 filler catalogs, one small table each: the scale a Databricks workspace showed
-- (the configurator used to pull every catalog's metadata before showing anything).
SET NOCOUNT ON;
DECLARE @i INT = 1, @name SYSNAME, @sql NVARCHAR(MAX);
WHILE @i <= 200
BEGIN
  SET @name = 'tenant_' + RIGHT('000' + CAST(@i AS VARCHAR(3)), 3);
  IF DB_ID(@name) IS NULL
  BEGIN
    SET @sql = N'CREATE DATABASE ' + QUOTENAME(@name) + N'; ALTER AUTHORIZATION ON DATABASE::' + QUOTENAME(@name)
      + N' TO ' + QUOTENAME(N'$(owner)') + N';';
    EXEC (@sql);
    SET @sql = N'USE ' + QUOTENAME(@name) + N'; CREATE TABLE dbo.events (event_id INT NOT NULL PRIMARY KEY, payload VARCHAR(64) NOT NULL); '
      + N'INSERT INTO dbo.events VALUES (1, ''' + @name + N' first event'');';
    EXEC (@sql);
  END
  SET @i += 1;
END
SELECT COUNT(*) AS catalogs_of_the_owner FROM sys.databases WHERE SUSER_SNAME(owner_sid) = N'$(owner)';
