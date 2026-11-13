# Connections — what the features leave out

## Still open

- core#1176 — Clone... saves the copy without the original's credentials (`cloneConnection` copies the parameters
  only); awaits a security ruling (a View-and-use grantee would get an editable copy of the owner's credentials).
  Claim when fixed: Clone... and OK on a fixture connection with a login made through the API, then the clone's
  login on the server — no outside database.

## Won't do

- TEST with a wrong host or credentials, special characters in Server or Db, port 0 or 65536: the database's or the
  connector's answer, not the dialog's.
- A name of spaces or with no letter in the New, Edit and Rename dialogs, and their CANCEL: the generic
  `Validators.notEmpty` and `DbEntityMeta.getRenameDialog` rules; CANCEL saving nothing is claimed by
  `connections-edit` (Clone...).
- A connection's queries deleted with it: ApiTests `Dapi: connection` 'deleting a connection deletes its queries'.
- Two connections of one name on a provider, and the Delete confirmation's text: as they are, by decision.
