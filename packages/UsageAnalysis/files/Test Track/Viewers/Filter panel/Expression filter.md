
1. Open demog
3. Open the **Filter Panel**
1. **Hamburger menu > Add filter > Expression filter**
4. Add filters:
   - **Numerical columns:**
     - `WEIGHT` `>` `50`
     - `HEIGHT` `<` `160`
   - **String columns:**
     - `SEX` `equals` `F`
     - `RACE` `contains` `an`
   - **Datetime columns:**
     - `STARTED` `after` `2005-01-01`
3. Switch between **AND/OR** mode
4. Right-click the filter value and select **Remove Query** for some filters
5. Click the **Free text** icon ( _I_ ).
2. Enter the following expressions: 
   * `AGE > 30 and SEX = M`
   * `AGE < 60 and HEIGHT > 190`
12. Save/apply the layout
---
{
"order": 5,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
