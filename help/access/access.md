<!-- TITLE: Access -->
<!-- SUBTITLE: -->

# Data access

Datagrok can access any machine-readable data source and retrieve either structured or unstructured data in a safe,
efficient, scalable, managed, and reproducible manner. Commonly supported data source types include
[databases](data-connection.md#connectors),
[file storages](file-shares.md), and web services. If anything is missing, it could be implemented as a platform
extension.

No matter what the data source type is, you would need to create a data connection in order to access the data. Each
data source type has its own connection parameters (for instance, that would be server/port/db for Postgres, and
region/bucket/folder for S3). Each connection might also be associated with login credentials (for instance, that would
be login and password for Postgres, and secret key for the S3 bucket). These credentials are encrypted and stored in a
special credential storage.

Connections are first-class entities in the Datagrok platform, and as such could be shared with users or groups (with
the specified privileges), included in projects, commented on, versioned, audited, and so on.

See also:

* [Data connection](data-connection.md)
