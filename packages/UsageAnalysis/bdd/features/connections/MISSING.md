# Connections edge cases — what the features cannot say yet

What the connection guards (folded into `connections-add`, `connections-edit` and `connections-delete`) leave out,
written for them: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one.

## Still open

### Decided and written

A clone of a connection carries its credentials, and a rename in the Edit dialog keeps them (decided 2026-10-27):
`connections-add` clones the stand's NorthwindTest and claims its schemas listed. Today the clone says
"Unavailable": core#1176, a known failure; the rename is claimed once the clone lists its schemas.

Claimed once core#1176 is fixed. The scenario ends at its pin, the Schemas node of the clone; with the tag
gone, its sweep becomes `Given no connection named "BDDConnClone{time}, BDDConnCloneRen{time}" is on the server`
and these steps go after the pinned one:

```gherkin
    Given Databases---Postgres---BDDConnClone{time}---Schemas tree node inside browse tree is expanded
    Then Databases---Postgres---BDDConnClone{time}---Schemas---public tree node inside browse tree should be visible
    When user picks "Edit..." from the context menu of Databases---Postgres---BDDConnClone{time} tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    When user enters "BDDConnCloneRen{time}" into Name input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 1 connection named "BDDConnCloneRen{time}" should be on the server
    Given Databases---Postgres---BDDConnCloneRen{time} tree node inside browse tree is expanded
    Then Databases---Postgres---BDDConnCloneRen{time}---Schemas tree node inside browse tree should become visible within 60 seconds
    Given Databases---Postgres---BDDConnCloneRen{time}---Schemas tree node inside browse tree is expanded
    Then Databases---Postgres---BDDConnCloneRen{time}---Schemas---public tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
```

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
