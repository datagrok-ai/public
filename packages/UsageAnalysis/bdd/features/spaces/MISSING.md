# Spaces edge cases — what the features cannot say yet

Written for the edge-case features `space-edge-cases` and `space-content-edge-cases`, with no library, binding or
code change. Each item below is a case those features leave out: the expected behaviour has no source yet, a
defect blocks the path, or the check belongs to another test layer.

## Still open

### Questions for owners

1. **A space deleted while it holds an entity the deleting user may not delete on its own (core#378).** A project
   that its owner may not delete was moved into a space the user may delete; deleting the space removed the
   project too. Should the delete be refused, keep the entity (moved back to its author's space), or delete it?
   The scenario is written once this is decided (the UI part: the refusal message, or the entity still listed).
2. **Rename a space to an empty name.** The Create Space dialog disables OK on an empty name (`spaces-create`);
   nothing says what the Rename dialog does. Expected: the same as create?
3. **A root space named like a system namespace ("Chem", "Bio") (core#358).** The ticket asks to allow it; today
   it is refused with "Root project with same name already exists". Not written yet: the feature's sweep deletes
   spaces by name, and a sweep of "Chem" must not reach the Chem package's namespace. Needs a decision on the
   expected name rules and a fixture name that is safe to sweep.
4. **A child space dragged onto the Spaces root.** `TestTrack/Spaces/spaces-missing-coverage.md` lists it as
   forbidden; datlas has "promote subspace to root", so it may be a legitimate promote. Which one?
5. **A package function dragged onto a space.** `spaces-missing-coverage.md` says it should be blocked;
   `help/datagrok/concepts/project/space.md` says entities from packages "can only be linked. Dragging them to a
   space automatically links them to it". The two sources disagree.

### Waiting for core#1042 (Link of a file into a space)

6. **The same file linked twice** (`TestTrack/Spaces/TESTS.md` #33): no error, and the space shows the file.
7. **A linked file gone with its source** (`TESTS.md` #35): the source file is deleted, and the link disappears
   from the target space. Copying a Demo file and then linking it from the source space is the old suite's path;
   both build on the Link that core#1042 breaks (`space-content-edge-cases`, the known failure).
8. **core#382 with a file**: a linked file whose original moved to another space. The project variant is claimed
   (green) in `space-content-edge-cases`.

### API or server tests to write

9. **Names that become storage paths and namespaces**: non-ASCII, "/", ":" and quotes in a space name, and files put
   into such a space stay readable (datlas `spaces/spaces_test.dart`, next to "create space with space in name and
   verify files").
10. **Two root names that differ only in case** ("Assays" and "assays"): the root check compares the name with `=`
    while `makeUniqueName` lower-cases; datlas, once the expected rule is stated.

## Won't do

- **core#692** (a View-and-use grantee links a query run into a shared space, then both users get an error on
  hover): since core#693 the Link is refused (`sharing-recipient`, "The Link is refused"), so the path that led to
  the error is closed.
- **Cascade of a delete to child spaces and copied files**: datlas `spaces_test.dart` ("delete space with child
  spaces cascade", "copy files between spaces"); the UI shows nothing more than the space gone.
