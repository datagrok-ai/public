1. Open demog
2. Open the **Filter Panel**
3. **Hamburger menu > Add filter > Expression filter**
4. For Expression filter:
   - Column set to `WEIGHT`, Operation set to `>`, in the Value input field enter 50, click the Add filter (+) icon
   - Column set to `HEIGHT`, Operation set to `<`, in the Value input field enter 160, click the Add filter (+) icon
   - Column set to `SEX`, Operation set to `equals`, in the Value input field enter F, click the Add filter (+) icon
   - Column set to `RACE`, Operation set to `contains`, in the Value input field enter an, click the Add filter (+) icon
   - Column set to `STARTED`, Operation set to `after`, in the Value input field enter 01/01/1991, click the Add filter (+) icon
5. On the Expression filter header, click OR icon - filters should use OR logic operation when applied
6. On the Expression filter header, click AND icon - filters should use AND logic operation when applied
7. Right-click the first created expression filter and select **Remove Query** - this filter should be removed
8. Hover over the Expression filter - icons appear on its header - click the **Free text** icon ( _I_ ) - the Search input field appears
9. In the Search input field, enter the `AGE > 30 and SEX = M` expression and click Enter
10. In the Search input field, enter the `AGE < 60 and HEIGHT > 190` expression and click Enter
11. Uncheck the first 4 expression filter rules - verify that filtered row count changes
12. Save the layout (Ctrl+S)
13. Close the Filter Panel
14. Apply the saved layout - the filter panel with unchecked rules should be restored
---
{
"order": 5,
"datasets": ["System:DemoFiles/demog.csv"]
}
