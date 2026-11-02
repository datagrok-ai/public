# Spaces edge cases — what the features cannot say yet

Written for the space edge cases, now @edge-cases scenarios of `spaces-rename`, `spaces-create`, `spaces-delete`,
`spaces-drag-and-drop` and `spaces-search`, their known failures in `spaces-known-failures`, with no library, binding or code change. Each item below is a case those scenarios leave out: the expected behaviour is not decided yet, a
defect blocks the path, or the check belongs to another test layer.

## Still open

### Waiting for a decision

1. **A space deleted while it holds an entity the deleting user may not delete on its own (core#378).** Decided:
   the delete of such a space is refused. How the refusal works is what core#378 is to design; the scenario (the
   refusal message, and the space and the entity still listed) is written once it is built.

### Waiting for core#1042 (Link of a file into a space)

**The linked file in the target space.** The known failure of core#1042 in `spaces-known-failures` ends at its pin, the
"Move entity" dialog closing. Claimed once core#1042 is fixed, after that step:

```gherkin
And no error or warning balloon should have been shown
When user double-clicks on BDDEdgeFileDst{time} tree node inside browse tree
Then the "BDDEdgeFileDst{time}" view should be current
And "edgelinkcars{time}.csv" link in gallery should be visible
```

2. **The same file linked twice into one space**: the second Link changes nothing, the space shows the file once
   and no error is logged (decided).
3. **A linked file gone with its source**: the source file deleted, the link goes from the target space too, a
   Link being the same entity (decided).
4. **core#382 with a file**: a linked file whose original moved to another space. The project variant is claimed
   (green) in `spaces-drag-and-drop`.

### Waiting for core#976 (a space deleted with what is linked into it)

**A package script dropped on a space is linked with no "Move entity" dialog** (`ProjectMeta.doDropEx`; a function of
a package, e.g. a BiostructureViewer script, can only be linked). No step removes a link from a space, and while
core#976 lasts the space's delete at feature end would delete the linked package script with it. Written once
core#976 is fixed (the space's cleanup then leaves the script), or once a step removes a link from a space.

### API or server tests to write

5. **Names that become storage paths and namespaces**: non-ASCII, "/", ":" and quotes in a space name, and files put
   into such a space stay readable (datlas `spaces/spaces_test.dart`, next to "create space with space in name and
   verify files").
6. **Two root names that differ only in case** ("Assays" and "assays"): the root check compares the name with `=`
   while `makeUniqueName` lower-cases; datlas, once the expected rule is stated.

## Won't do

- **core#692** (a View-and-use grantee links a query run into a shared space, then both users get an error on
  hover): since core#693 the Link is refused (`sharing-recipient`, "The Link is refused"), so the path that led to
  the error is closed.
- **Cascade of a delete to child spaces and copied files**: datlas `spaces_test.dart` ("delete space with child
  spaces cascade", "copy files between spaces"); the UI shows nothing more than the space gone.
- **Escape on the space dialogs**: every dialog closes on Escape as on CANCEL (`Dialog._closeOnEsc`), and CANCEL is
  claimed by spaces-create, spaces-rename, spaces-delete and spaces-drag-and-drop. The one branch of its own, Escape
  closing a popup over a dialog first, is claimed in `sharing-dialog`.
