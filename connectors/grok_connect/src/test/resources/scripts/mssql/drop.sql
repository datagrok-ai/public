EXEC sp_msforeachtable "ALTER TABLE ? NOCHECK CONSTRAINT all";
EXEC sp_MSforeachtable @command1 = "DROP TABLE ?";