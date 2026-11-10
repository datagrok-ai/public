# Spaces — what the features do not claim yet

## Still open

- core#976 — deleting a space deletes a project that was only linked into it; claim when fixed: a datlas
  `spaces_test.dart` case (link a project, delete the space, the project is still found with its namespace).
- core#378 — a space holding an entity its deleter may not delete is to refuse the delete (decided); the refusal is
  core#378's design, claimed once it is built.
- A package script dropped on a space is linked with no "Move entity" dialog (`ProjectMeta.doDropEx`): waits for
  core#976 (the space's cleanup would delete the script) or for a step that removes a link from a space.

## Won't do

- core#692 (a View-and-use grantee's link of a query run into a shared space): the Link is refused since core#693
  (`sharing-recipient`, "The Link is refused").
- How far a delete reaches (child spaces, copied files): datlas `spaces_test.dart`; the UI shows only the space gone.
- Escape on the space dialogs: every dialog closes on Escape as on CANCEL (`Dialog._closeOnEsc`).
- A child moved out to the root onto a taken root name: datlas 'forbid duplicate root space names';
  `spaces-create` claims the refusal balloon on Create Space.
- A link that follows its original moved to a third space (core#382) on spaces of its own: claimed on a project in
  `projects-move` (BDDMoveStay's link in the second space); not on a file.
- Linking edge cases (the same file twice, a link whose source is deleted): the server's link semantics.
- Characters of the search syntax typed into the spaces search: the filter-grammar corpus
  (`grok_shared/test/filter_grammar_corpus_test.dart`).
- A nested favorite also listed at the top level of Favorites: a product choice (7be502d464), not guarded.
- Space names as storage paths beyond spaces (non-ASCII, "/", ":", two roots differing in case): datlas
  'create space with space in name and verify files' covers the path rule.
