# Apps edge cases — what the features cannot say yet

Written for `app-address-edge-cases.feature`, with no library, binding or code change. The feature claims that an
address naming no app says "Application not found" (`routing.dart` `_appsHandler`). What follows has no source
for its expected behaviour.

## Still open

### Questions for owners

1. **The view after "Application not found".** On localhost (2026-10-26) the address `/apps/<no such app>` opened
   the Apps view in front of the table that was current, after the error balloon. Is the Apps view the intended
   landing place, or should the current view stay?
2. **An existing app's address with a path that names nothing inside it** (`/apps/Tutorials/Demo/<no such demo>`):
   the path goes to the app's own routing, and no source says what the app should show.
3. **Hit Design: Add new row, then CANCEL in the molecule sketcher.** Is the empty row kept or removed, and is a
   V-iD registered for it? The campaign is shared, so the answer decides whether a feature can run it at all.

## Won't do

- **Edge inputs in the demos and guides.** They demonstrate a happy path (`guides/INDEX.md`); edge data in a demo
  would test the viewers, which have their own edge-case features.
