# PeptiHit

PeptiHit is Hit Design for peptides. A design is a sequence in HELM notation: you enter sequences in a
spreadsheet, PeptiHit converts each one to its molecule, calculates the molecule's properties, gives it a V-iD,
and keeps the designs of a campaign in stages until you submit the ones you pick. Campaigns are saved and
shared with others, as in Hit Design.

PeptiHit is in `Browse` → `Apps` → `Peptides`. It needs the Bio and Helm packages: Bio converts sequences to
molecules, Helm draws and edits the HELM cells.

The application is built around templates and campaigns. A template says what every campaign made from it
records and calculates; a campaign is the table of designs itself.

## Templates

A PeptiHit template holds the same settings as a Hit Design template, and is created the same way, with
the `+` (New Template) or the copy icon (Clone template) beside the `Template` choice on the landing:

- **Name** : Identifies the template.

- **Key** : The prefix of the campaign codes (`PHD-1`, `PHD-2`, ...).

- **Layout** : The layout the campaign table opens with. PeptiHit applies the template's layout to every
  campaign that has no layout of its own, and `Apply template layout` on the ribbon applies it again at any time.
  A layout picked from a mounted file share is read again from that file each time it is applied.

- **Additional fields** : Fields asked for when a campaign starts (a chemist's name, a target, a deadline),
  each with a type and whether it is required.

- **Stages** : The lanes of the Progress Tracker. A new row starts in the first stage.

- **Compute** : The calculations every design gets: molecular descriptors and the functions tagged
  `HitTriageFunction` in the installed packages. They run on the molecule a sequence converts to, not on the
  sequence.

- **Submit function** : The function the campaign's rows are submitted to (tagged `HitTriageSubmitFunction`).

The package ships a demo template, `Peptide hits Demo` (key `PHD`), with a `Chemist` field, four stages
(Design, Consideration, Synthesis, Testing) and three descriptors (MolWt, HeavyAtomMolWt, FractionCSP3).

## Campaigns

- **Creation** : Pick a template, fill in its fields and an optional campaign name, and press `Start`.
  The campaign opens with one empty row and the HELM editor on it.

- **Designs** : The `Helm` column holds the sequences. When a sequence is entered or changed, PeptiHit converts
  it to its molecule (Bio's To Atomic Level, nonlinear sequences included) into the `Molecule` column, runs the
  template's compute on that molecule, and registers the molecule's V-iD. V-iDs come from one dictionary shared
  with Hit Design: one structure has one V-iD.

- **Editing** : Add rows with the `+` icon on the ribbon or `Add new row` in a sequence cell's context menu;
  `Duplicate Sequence` copies a row's sequence into a new row, and `Re-Run Calculations` converts and calculates
  a row again.
  `Calculate additional properties` adds calculations to the campaign and runs them over the whole table.

- **Saving** : The campaign is saved as you work: after each calculation, stage change, merge, permissions
  change and Save in the Submit dialog; each save keeps the table's layout too. Two people can work on one
  campaign at once: a save merges in what the other session saved meanwhile, and the `Refresh` icon tells you
  when there is something new to read.

- **Collaboration** : A campaign's address (`?campaignId=PHD-1`) opens it directly, for anyone who may view it.

PeptiHit has no designer function button: the `HitDesignerFunction` functions are offered in Hit Design only.

## Getting started

The landing lists the campaigns you may view under `Continue Campaign`: open one by its name. The table can be
grouped by template, status, author, last modifier or a campaign field, its columns chosen with the eye icon,
and refreshed from the server with the sync icon. A campaign you may only view opens with "View only — changes
are not saved" in its status bar.

Start a new campaign under `New Campaign`: pick the template, fill in its fields, and press `Start`.

In the campaign, type or draw sequences in the `Helm` column. The `Progress Tracker` button on the ribbon opens
the stages as lanes: drag a card to another lane to move the design to that stage; the moves are saved when the
tracker closes. `Modify Stages` in the tracker adds, removes, renames and reorders stages; rows in a stage that is
removed move to the first one.

`Merge another table into the campaign` brings rows and columns of another table (a file from your computer or a
file share) into the campaign, matched by canonical SMILES or by V-iD.

`Submit...` on the ribbon shows what the campaign is and lets you change its name, its folder and its status (you
are offered the statuses other campaigns use). `Save` stores the status; `Submit` runs the submit function on
the campaign's rows (the filtered ones) and then saves.

## Adding custom compute and submit functions

PeptiHit uses the same compute and submit functions as Hit Design, written in any installed package. A compute
function takes the table and its molecule column first; in PeptiHit that column is the `Molecule` column the
sequences convert to. A submit function takes the table and the name of the molecule column. See
[Hit Design](README_HD.md#adding-custom-compute-and-submit-functions) for the signatures and for examples in
TypeScript, JavaScript, Python and SQL.

## Permission Management

Each campaign keeps who may view it and who may edit it, by user or group. A new campaign is shared with the
groups in the HitTriage package settings (all users by default).

To change a campaign's access, use `Manage campaign permissions` in its row of the campaigns table, or the users
icon on the campaign's ribbon. Can Edit includes Can View. With no one left under Can Edit, only the campaign's
author can change it; with no one under Can View, all users can view it.

To change the default access of new campaigns:

1. Go to `Browse` → `Platform` and click on `Plugins` (or directly via a link at `{your_datagrok_url}/packages`).
2. Under Cheminformatics, click on `HitTriage` package.
3. On the context panel, expand the `Settings` pane.
4. Modify the `Default view permissions` and `Default edit permissions` fields.

![peptiHitReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/edit-default-hd-sharing.gif?raw=true)
