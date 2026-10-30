# Bio transform features — what is still missing

## Bio panes

Feature: `get-region-pane.feature` (30 Oct 2026, local stand).

### No such pane (Bio panes)

- **An MSA context pane.** `Bio:multipleSequenceAlignmentDialog` is declared with `@grok.decorators.panel({name:
  'MSA', tags: ['bio', 'panel'], 'top-menu': 'Bio | Analyze | MSA...'})`, but it takes no input, so the platform
  never offers it as a pane: on a fasta column made current, the context panel shows Details, Filter, Actions,
  Colors, Style, Settings, Plots, Advanced, Bioinformatics (Get Region, Manage Monomer Libraries, Sequence
  Renderer) and Peptides, and the Actions pane does not list MSA either (measured 30 Oct 2026 on the local stand).
  The MSA dialog is `analyze/msa.feature`. Is the `panel` decorator on the dialog intended (a pane was planned),
  or left over? Nothing to write until a pane exists.
