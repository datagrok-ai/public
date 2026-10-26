# Connections edge cases — what the features cannot say yet

What `connection-edge-cases.feature` leaves out, written for it: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one.

## Still open

### Questions for owners

1. **Edit... that changes only the name keeps the credentials.** The password field comes back masked, and
   `files/TestTrack/Connections/COVERAGE.md` notes a credential-cache lag in edit. Is the contract that
   credentials not touched in the dialog stay as they were (a query of the connection still runs)?
2. **Two connections with one name on one provider.** The server suffixes the grok name (`_1`) but keeps the
   friendly name, so the tree shows two identical nodes (`makeUniqueName`, `repository_query.dart`). Refused,
   renamed, or allowed? Once decided, an ApiTests case (the class of core#546).
3. **Port 0, 65536 or letters in the New connection dialog.** The only numeric field; its validation is not
   stated anywhere.
4. **Clone... then OK.** Only Clone + CANCEL is claimed (`connections-edit.feature`). Does the copy carry the
   credentials?
5. **The Delete confirmation of a connection that has queries.** It reads `Delete connection "<name>"?` and the
   queries go with it (help `access/databases/databases.md`, the caution under Query Editor, claimed in
   `connection-edge-cases.feature`). Should it name the queries it will delete?

### API or server tests to write

6. **A connection delete removes its queries.** No ApiTests case checks it (`dapi/connection.ts` deletes
   connections without queries). The feature claims it through the UI delete; the server rule belongs in
   ApiTests `dapi/connection.ts`.
7. **core#607, credentials row and wrapper project left after a delete** — already requested in
   `projects/MISSING.md` §3 (datlas).

## Won't do

- **TEST with a wrong host or wrong credentials.** What TEST answers is the database server's reply to a login,
  not the dialog's (the bdd library's CLAUDE.md, "What never becomes a feature"; the feature headers of
  `connections-add` and `connections-edit` say the same).
- **Special characters in Server or Db.** A connector concern, for grok-connect's tests.
