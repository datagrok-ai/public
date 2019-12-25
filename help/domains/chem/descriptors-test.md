<!-- TITLE: Tests: Chem Descriptors -->
<!-- SUBTITLE: -->

# Tests: Chem Descriptors

[Molecular descriptor](../domains/chem/descriptors.md) is the final result of a logic and mathematical procedure which
transforms chemical information encoded within a symbolic representation of a molecule
into a useful number or the result of some standardized experiment.

## Testing scenario

1. Open *smiles.csv* file

1. Click on *"smiles"* column

1. Open *"Descriptors"* dialog in *"Active"* tab on [Property Panel](../features/property-panel.md) 
   (or use context menu for *"smiles"* column, **Chem | Descriptors...**)
   * *"Descriptors"* dialog is open

1. Mark "EState VSA" group in the tree
   * "EState VSA" group is marked in the tree
   * "EState VSA" group is in bold
   * Next to the group name is number of marked [descriptors](../domains/chem/descriptors.md) (21/21)

1. Expand "EState VSA" group

1. Deselect "EState_VSA1" [descriptor](../domains/chem/descriptors.md)
   * "EState VSA" group is no longer shown marked in the tree
   * "EState VSA" group remains bold
   * Number of selected indicators is shown as (20/21)

1. Expand "Lipinski" group

1. Mark [descriptors](../domains/chem/descriptors.md) *"FractionCSP3"*, *"HeavyAtomCount"* and *"NHOHCount"* in "Lipinski" group
   * When [descriptor](../domains/chem/descriptors.md) is marked, "Lipinski" group is allocated with bold
   * Number of selected [descriptors](../domains/chem/descriptors.md) changes next to group name

1. Execute dialog with previously specified parameters
   * Status bar with progress is displayed
   * Columns with selected [descriptors](../domains/chem/descriptors.md) added to table for each row

1. Open *"Descriptors"* dialog for *"smiles"* column

1. Select last entry from history
   * 20 [descriptors](../domains/chem/descriptors.md) in "EState VSA" group and 3 [descriptors](../domains/chem/descriptors.md) of  "Lipinski" group are marked
   * "EState VSA" and "Lipinski" groups are  bold and counter of selected [descriptors](../domains/chem/descriptors.md) is shown for each

1. Test non-functional modules (UI, help, navigation, console)

See also:
 * [Molecular descriptor](../domains/chem/descriptors.md)
 * [Substructure Search Test](../tests/substructure-search-test.md)
 * [Get R Groups Test](../tests/get-r-groups-test.md)
