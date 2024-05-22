Data: demog.csv
1. Add new columns with the following formulas:

* if(true, "yes", error("Error"))
* if(length(${DEMOG}) > 1,substring(${DEMOG}, 0, 3), "no")

(TODO: Add expected result. Should the new column have some exact data?)

---
{
  "order": 4
}
