# Campaign catalog

A u2 `CatalogTable` of the campaigns, with the columns Name, Code, Kickoff and Owner. Kickoff reads the way dates
read across the platform, year/month/day without padding (`2026/6/29`); Owner reads the owner's friendly name. The
user can group the list by Kickoff or by Owner, through groupings named `Kickoff` and `Owner`.

The plugin calls:

```ts
type Campaign = {code: string, name: string, kickoff: Date, owner: {login: string, friendlyName: string}};
export function campaignCatalog(campaigns: Campaign[]): CatalogTable<Campaign>
```
