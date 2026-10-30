# Filter Builder feature — what is still missing

## Compare Columns, Filter Builder

Feature: `filter-builder.feature` (30 Oct 2026, local stand).

### Questions (Filter Builder)

- **No help page** covers the Filter Builder card in the filter panel (`help/` has no "filter builder"); the
  feature claims what `filter-builder-filter.ts` does. The U2Demo catalog features cover the u2 control itself,
  not this host.
- **The builder/query switch is an icon with a title only.** It is reached as `"Switch between the builder and the
  query text" icon` (its title); it carries `data-u2-part="query-toggle"` outside the filter builder's root, so
  no part phrase of the card reaches it. Should it become a part of the card (or take an aria-label)?
- **The card's own on/off checkbox and its x** are the filter panel's, shared with every card; not claimed here.
