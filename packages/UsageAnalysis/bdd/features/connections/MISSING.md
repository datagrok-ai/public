# Connections edge cases — what the features cannot say yet

What `connection-edge-cases.feature` leaves out, written for it: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one.

## Still open

### Decided and written

A clone of a connection carries its credentials, and a rename in the Edit dialog keeps them (decided 2026-10-27):
`connection-edge-cases` clones the stand's NorthwindTest and lists its schemas before and after the rename. Today
the clone says "Unavailable": a known failure with no ticket yet.

### API or server tests to write

1. **A connection delete removes its queries.** No ApiTests case checks it (`dapi/connection.ts` deletes
   connections without queries). The feature claims it through the UI delete; the server rule belongs in
   ApiTests `dapi/connection.ts`.
2. **core#607, credentials row and wrapper project left after a delete** — already requested in
   `projects/MISSING.md` §3 (datlas).

## Won't do

- **TEST with a wrong host or wrong credentials.** What TEST answers is the database server's reply to a login,
  not the dialog's (the bdd library's CLAUDE.md, "What never becomes a feature"; the feature headers of
  `connections-add` and `connections-edit` say the same).
- **Special characters in Server or Db.** A connector concern, for grok-connect's tests.
- **Two connections with one name on one provider**: kept as it is, by decision (the grok name gets `_1`, the
  friendly name stays).
- **Port 0, 65536 or letters in the New connection dialog**: the port may also come in the server name or the
  connection string, so the field is not bounded; TEST is the check, by decision.
- **The Delete confirmation naming the queries it deletes**: the present text is enough, by decision.
