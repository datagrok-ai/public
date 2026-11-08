# Sharing: what the features do not claim yet

Who may see, change, re-share or lose a shared entity is the server's answer: ApiTests `Dapi: permissions`
(`src/dapi/permissions.ts`), not a feature.

## Still open

- An entity's own project (`projects.is_entity`) and the grants on it outlive the deleted query, script or connection
  (`DELETE /api/projects/<id>` keeps them): needs a product fix before the features can leave the stand clean.
- The Share dialog's "This entity will also be shared to view and use" notice: needs a fixture whose project holds a
  second entity (a query saved through the UI, a script in a project).
- A step setting the level of an existing grant row (`user sets the access level of the sharing user in {element}
  to {string}`); only the new grant's selector is reachable today.
- Chats pane readings: a message's author, the order of messages, a message's own menu, a mention.
- ApiTests `Dapi: permissions`: a folder share gives read but not write, and files outside the folder stay closed.
- A disabled user is not offered in the suggestions: needs a disabled-account fixture step and a suggestions reading.

## Won't do

- A grant revoked while the recipient has the entity open: two sessions at once, which one page cannot hold.
- What a View-and-use recipient may do beyond re-sharing: ApiTests `permissions.ts`; `sharing-recipient` claims the UI.
- Escape on the Share dialog: it closes as CANCEL does (`Dialog._closeOnEsc`).
- A grant removed in the dialog and then cancelled: every dialog writes on OK only; the removal itself is claimed in
  `sharing-recipient` and EDA `share-model`.
