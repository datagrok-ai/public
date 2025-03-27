---
title: "Logo Summary Table"
---

The Logo Summary Table shows WebLogo, activity distribution and statistics for each sequence cluster in your data. Statistics include mean difference, p-value, count and ratio. You can add other column aggregations to the table from its property panel. Numeric columns will use your chosen numeric aggregation (average by default) and categorical columns will add a pie chart showing the distribution within a cluster.

## Customization

You can modify the Logo Summary Table by changing the following properties from the property panel:

- **Sequence**: The column containing the sequence data.
- **Clusters**: The column containing the cluster data.
- **Activity**: The column containing the activity data.
- **Activity Scaling**: The scaling method for the activity data.
- **WebLogo Mode**: The mode for the WebLogo, Entropy or 100%.
- **Members Ratio Threshold**: Internal ratio threshold to filter out smaller clusters.
- **Aggregations**: Additional columns to aggregate in the table.

## Responsiveness

The Logo Summary Table reacts to filters applied to the data and shows only clusters of filtered rows. When filters are active, the WebLogo, activity distribution and statistics are calculated using only the filtered rows.

## Interactivity

Selecting a cluster will select every sequence in that cluster.

The Logo Summary Table follows Datagrok's standard selection pattern:

|                  |                               |
|------------------|-------------------------------|
| Click            | Single selection of Cluster   |
| Ctrl+Click       | Toggle Cluster selection state|
| Shift+Click      | Add Cluster to selection      |
| Ctrl+Shift+Click | Deselect Cluster              |

Hovering over cluster cell will show tooltip with distribution and statistics of the cluster and highlight corresponding
points on other viewers.

[Logo Summary Table](./img/LST.gif)
