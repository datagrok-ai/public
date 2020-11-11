<!-- TITLE: Shape Map -->
<!-- SUBTITLE: -->

# Shape Map

Shows a map that is applicable for the specified dataset. Typically, it would 
represent a geographical area (countries, states, counties, etc), but it also
supports arbitrary shapes (such as a store floor plan, brain regions,
or EEG electrodes).

When opened, a viewer automatically determines the best map that is applicable 
to the current dataset. If more than one map fits the data, it is chosen
arbitrarily; right-click and choose the one you want from the pop-up menu.

Geographical regions can be represented by different names (such as Pennsylvania and PA), 
and there is an elaborate system in place that understands synonyms and abbreviations. However,
sometimes it can't map names to the regions. To identify these records, select
"Select not matching rows" from the popup menu.

![Shape Map](../../uploads/viewers/shape-map-pa-counties.png "Shape Map")
![Shape Map](../../uploads/viewers/shape-map-plate.png "Shape Map")

## Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/7MBXWzdC0-I?start=3650" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also:
* Applicable tables: #{x.demo:germany_grp_by_state}, #{x.demo:pa_income_by_county}, #{x.demo:ua_population}  
* [Viewers](../viewers.md)
* [JS API: Shape Map](https://public.datagrok.ai/js/samples/ui/viewers/shape-map)
