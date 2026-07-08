# Peptide Generation

The Peptide Generation viewer designs new candidate peptides from your data and ranks them by predicted activity. It reuses the same per-monomer-per-position statistics as the [Sequence Variability Map](sequence-variability-map.md) and [Most Potent Residues](most-potent-residues.md) viewers, so it fits naturally into a [Peptides SAR](../../datagrok/solutions/domains/bio/peptides-sar.md) analysis. It is an optional viewer and is not added by default. Add it like any other viewer from the **Add viewer** menu.

Candidates are scored with an additive positional (Free-Wilson) model:

```
predicted activity = global mean + Σ over positions (mean activity of the chosen monomer at that position − global mean)
```

Each per-position contribution is damped by an empirical-Bayes shrinkage factor so that poorly-supported monomers cannot inflate the prediction. A bounded beam search combines the best-scoring monomers at each position into whole peptides toward your chosen High or Low activity target, and the viewer keeps the top peptides.

:::note

Predicted activity is a ranking score and an optimistic estimate, not a calibrated measurement. Because contributions add across positions, it can extrapolate beyond the observed activity range, and the model treats positions as independent (it cannot see epistasis). Use the Confidence, Min support, Significant positions, and basis columns to judge how much to trust each row.

:::

## Result columns

Each generated peptide is one row, sorted by predicted activity:

- **Peptide**: The generated sequence, rendered as a macromolecule.
- **Predicted activity**: The additive-model estimate, used as the sort key.
- **Confidence**: The fraction of varied positions whose chosen monomer is statistically significant.
- **Min support**: The smallest number of observed sequences behind any chosen monomer.
- **Mean p-value**: The average p-value across the chosen monomers.
- **Significant positions**: The number of positions whose chosen monomer is significant.
- **Novel**: Whether the peptide is absent from the source data.

## Basis columns

The trailing columns, one per position, show the basis for each prediction. Every cell displays the chosen monomer, an effect circle, and the signed contribution:

- **Monomer**: Colored with its [monomer library](../../datagrok/solutions/domains/bio/bio.md#manage-monomers) color and fitted to the column width.
- **Circle size**: The effect magnitude (absolute contribution), normalized within the row. A bigger circle means a stronger effect on activity.
- **Circle color**: The direction and significance of the effect. Red raises activity, blue lowers it, a pale color means weak significance, and grey means there was too little data to run the test.
- **Number**: The signed contribution to the predicted activity.

Hovering over a basis cell shows the activity distribution of sequences with that monomer at that position versus the rest, together with the statistics behind the choice.

## Customization

You can modify the Peptide Generation viewer by changing the following properties from the property panel:

- **Sequence**: The column containing the sequence data.
- **Activity**: The column containing the activity data.
- **Activity scaling**: The scaling method for the activity data.
- **Activity target**: Whether to design peptides that maximize (High) or minimize (Low) activity.
- **Peptide count**: How many peptides to generate.
- **Candidates per position**: The number of best-performing monomers considered at each position.
- **Beam width**: The beam search width. Larger values explore more combinations at higher cost.
- **Min support**: The minimum number of observed sequences backing a monomer for it to be trusted.
- **Max p value**: Only consider monomers with a p-value at or below this value (1 disables the filter).
- **Shrinkage strength**: Damps poorly-supported monomers by a factor of count / (count + k). Higher values are more conservative, 0 disables shrinkage.
- **Exclude existing**: Exclude generated peptides that already exist in the source data.
- **Only significant positions**: Vary only positions that have at least one significant monomer, and keep the most common monomer elsewhere.

## Responsiveness

The Peptide Generation viewer reacts to filters applied to the data. When filters are active, the statistics, predictions, and generated peptides are computed using only the filtered rows. When available, it shares the cached statistics of the other SAR viewers to avoid recomputation.

## Interactivity

Hovering over a basis cell displays a tooltip with the activity distribution and statistics of the corresponding monomer-position and highlights the matching rows in other active viewers.

[Peptide Generation](../../datagrok/solutions/domains/bio/img/peptides/peptide-generation.png)
