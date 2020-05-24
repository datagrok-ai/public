<!-- TITLE: View Layout -->
<!-- SUBTITLE: -->

# View Layout

View Layout contains relative positions of [viewers](../visualize/viewers.md) in a [table view](../overview/table-view.md),
along with the viewers' properties. By separating layouts from the actual data displayed, we now can
save current layout (**View | Layout | Save to Gallery**) and later apply it to a different dataset
(**View | Layout | Open Gallery**). 

To clone current view, either do **View | Layout | Clone**, or click on the plus sign on the view header strip, 
and choose **Clone**.

# Viewer Layout

Similarly to view layout, it is possible to reuse the settings of individual [viewers](../visualize/viewers.md).   
Under the 'Viewer' menu (either right-click on the viewer, or click on the 'hamburger' menu on the top left),
you will find the following layout-specific commands:
* Clone
* Save to Gallery

To open a gallery containing viewer and view layout suggestions: **View | Layout | Open Gallery**

## Layout Suggestions

Moreover, the platform will proactively suggest layouts (even created by
other users) to be applied to your data, if a previously saved layout can be used for visualizing currently
opened dataset. The suggestions are based on the [layout applicability](#layout-applicability), 
popularity and specificity of the layout, as well as on the previously observed actions.
 
To check layouts applicable to the current table, open 'Layouts' pane on the left. Alternatively,
open the specialized view/viewer layout suggestion panel: **View | Layout | Open Gallery** 

## Layout Applicability

When a view layout is saved, the visual arrangement of the viewers along with the metadata of the 
layout data columns (columns selected on viewers, such as "X" column on a scatter plot) gets saved. 
Metadata includes column name, type, semantic type, and other metadata.

To determine whether a view layout is applicable to a particular table, the platform checks whether
all layout data columns could be mapped to the table columns. There are few different ways a 
layout column gets mapped to the table column, which are evaluated in the following order:

1. Column names AND column types match
2. Both columns have the same [layout-id](../discover/tags.md#layout-id) 
3. Both columns have the same [semantic-type](../discover/tags.md#semantic-type) 
 
Columns get mapped in the following way
in this order:
1) 


See also:
* [Table view](../overview/table-view.md)
* [Self-learning platform](../learn/self-learning-platform.md)
