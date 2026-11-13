# Browse — what the features leave out

## Still open

- core#1152 — a `/files/` address naming no folder opens an empty folder view (logs "File or Directory … not
  exists"), a `/script/` address naming no script an empty Script view (logs the 404); the outcome is undecided, and
  the router's NullError on the next address followed is a bug either way. Claim when fixed: both addresses followed
  in `browse-routing`, what the fix shows, and the next address opening.

## Won't do

- Addresses that name nothing: `/p/` and `/apps/` are claimed in `browse-routing`, `/func/` in xamgle
  `routing_test.dart`; a deleted project's or space's address is the same lookup as a never-saved one, and `/s/`
  just lists the spaces.
- Files Browse cannot preview (an unknown extension, an empty .csv or .txt): Download is claimed on a CSV; an empty
  file's preview is the reader's.
- A file over the 10 MB preview limit, an empty folder, an entity deleted by another session: a large fixture on
  every run, nothing to claim, two sessions at once.
- core#215 (an entity readable inside a parent that is not): a permission matrix, ApiTests `Dapi: permissions`.
- A package's folder viewer (Clinical Case's SDTM launcher): the package is not on the BDD stand, and a folder viewer
  loads its package on every folder opened in Files.
