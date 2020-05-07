<!-- TITLE: Tests: Cluster Data -->
<!-- SUBTITLE: -->

# Tests: Cluster Data

[Clustering](cluster-data.md) is the task of grouping a set of objects in such a way that objects in the same group
(called a cluster) are more similar (in some sense or another) to each other than to those in other


## Testing scenarios:

1. Open "demog" dataset
  
1. Open "Cluster..." from  **Tools | Data Science** 

1. Add "Age", "Height" and "Weight" columns to field "Features"

1. Enter the number of [clusters](cluster-data.md) equal to 3 in the field "Clusters"

1. Set the value of field "Show scatter plot" as true
   * Column "clusters" was added to "demog" table
   * Viewer "[Scatter Plot](../visualize/viewers/scatter-plot.md)" was created where [clusters](cluster-data.md) are marked with color

1. Click on "Cancel" button in "Cluster" dialog
   * Column "clusters" and [scatter plot](../visualize/viewers/scatter-plot) disappeared
   
1. Add  all available columns to field "Features"

1. Click "OK" button
   *   Column "clusters" was added to "demog" table 
     
1. Test non-functional modules (UI, popup menu, help, navigation, properties, etc.)
   * Non-functional modules work correctly and are intuitive

   
See also:
  * [Cluster data](cluster-data.md)
  * [Scatter plot](../visualize/viewers/scatter-plot.md)
  * [Cluster](cluster-test.side)