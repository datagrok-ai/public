<!-- TITLE: Anonymize data -->
<!-- SUBTITLE: -->

# Anonymize data

Sometimes, a need arises to prepare a dataset that conveys an idea of the structure in the data along with the patterns
in it, but does not contain the real data points. This is what 'Anonymize Data' functionality is for.

The data in the selected columns gets replaced with the anonymized values. The new value depends on the column type, as
well as on the old value:

* **Categorical.** If the categories are defined for the specified column, the value is replaced with a random category.
  Otherwise, the value is replaced with the synthetic category (such as '
  race 5' for the 'race' column). To choose a category, click on the corresponding cell in the '
  category' column. Note that a table containing exactly one string column (representing a category)
  has to be imported before opening data anonymization dialog.
* **Numerical.** The number gets randomly changed. The scale of the change depends on the specified _number
  randomization factor_ **r**. The actual formula is:

  ```
  newValue = col[i] * math.pow(10, 2 * r * rnd.nextDouble() - r)
  ```

For selecting a random subset of the (possibly already anonymized) dataset,
see [Selecting random rows](../explore/select-random-rows.md).
