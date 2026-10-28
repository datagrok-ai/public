# Apps edge cases — what the features cannot say yet

Written for `app-address-edge-cases.feature`, with no library, binding or code change.

## Decided and written

The answers of 2026-10-27 are scenarios of `app-address-edge-cases`: after "Application not found" the Apps view is
current. A path naming no demo inside the Demo app is left out (2026-10-28): the Demo app opens after a while on
dev, while a local stand kept an empty view for two minutes.

## Won't do

- **Edge inputs in the demos and guides.** They demonstrate a happy path (`guides/INDEX.md`); edge data in a demo
  would test the viewers, which have their own edge-case features.
- **What happens inside an app** (Hit Design's rows, an app's own dialogs): each app's behaviour is its own task,
  outside this round.
