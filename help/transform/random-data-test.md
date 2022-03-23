<!-- TITLE: Tests: Random data -->
<!-- SUBTITLE: -->

# Tests: Random data

Adds a numerical column with [random data](random-data.md) with the specified distribution. Parameters of the
distribution can be edited as well.

## Testing scenarios

1. Open "demog" dataset

1. Open "Random data" dialog from  **Tools | Data Science**

1. Leave the values of all fields by default. Set the value of field "Show histogram" as true

* Column "Normal" was added to "demog" table
* Viewer "[Histogram](../visualize/viewers/histogram.md)" with data from "Normal" column was created
* Distribution of "Normal" column values corresponds to normal (you can evaluate it on a histogram)

1. Set values of fields "Seed" = 100, "Mean" = 10 and "Std. = 5" in "[Random data](random-data.md)"
   dialog

* Values in "Normal" column was changed corresponds to new parameters
* [Histogram](../visualize/viewers/histogram.md) was changed corresponds to new values of "Normal"
  column

1. Click on "Cancel" button in "Random data" dialog

* Column "Normal" and [histogram](../visualize/viewers/histogram.md) disappeared

1. Repeat steps 2-5 for all available distributions. For step #4, use the valid values of the available fields for each
   distribution

1. Open "Random data" dialog. Leave the values of all fields by default. Set the value of field "
   Show histogram" as true. Alternately select all the possible values ​​of the "Distribution"
   field ("Normal", "Log-Normal", "Binomial", "Poisson", "Uniform")

* When you select a new value for the "Distribution" field, a new column is added to "demog" table
* Names and values ​​of new columns correspond to the selected distribution
* [Histogram](../visualize/viewers/histogram.md) is created for each new column
* When "Distribution" field value changed, then parameters that are not suitable for the corresponding distribution
  become unavailable for editing

1. Click "OK" button

* New columns added to the "demog" table
* [Histogram](../visualize/viewers/histogram.md) is added to view for each new column

1. Test non-functional modules (UI, popup menu, help, navigation, properties, etc.)

* Non-functional modules work correctly and are intuitive

See also:

* [Random data](random-data.md)
* [Histogram](../visualize/viewers/histogram.md)
* [Random data Auto Test](random-data-test.side)
