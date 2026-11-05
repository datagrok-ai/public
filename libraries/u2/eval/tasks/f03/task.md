# Campaigns home

The app's home lists the user's campaigns in a u2 `CatalogTable`, with the columns Name, Code, Status, Kickoff and
Owner. `load()` fetches them, which takes a moment, so the home shows that it is loading meanwhile.

- Kickoff reads the way dates read across the platform, year/month/day without padding (`2026/6/29`); Owner reads
  the owner's friendly name.
- The list starts grouped by Status, each group headed by its status and how many campaigns it holds. The user
  can group it by Kickoff or by Owner instead, through groupings named `Status`, `Kickoff` and `Owner`.
- As in the platform's other lists, the row under the pointer is highlighted.
- A click on a campaign's name shows its details beside the list: **Code**, **Name**, **Owner** and **Structure**,
  captions on the left and values on the right as in a form. The structure is the platform's drawing of the
  campaign's reference compound, `grok.chem.drawMolecule(campaign.reference.smiles, 160, 100)`, an element of
  160×100 pixels. Under them, **OPEN** calls `open(campaign)`.

The plugin calls:

```ts
type Campaign = {code: string, name: string, status: string, kickoff: Date,
  owner: {login: string, friendlyName: string}, reference: {id: string, smiles: string}};

export function campaignsHome(load: () => Promise<Campaign[]>, open: (campaign: Campaign) => unknown): Control
```
