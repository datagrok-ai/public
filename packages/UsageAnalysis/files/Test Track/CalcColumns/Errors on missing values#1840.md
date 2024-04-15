Data: ForDGWithSigns.csv
1. Add new columns with the following formulas:

* Qnum((100-${withSigns2})/100, if(qualifier(${withSigns2})==">", "<", "="))
* Qnum((100-${withSigns2})/100, if(isnotempty("${withSigns2}"), if(qualifier(${withSigns2})==">", "<", "="), null))
---
{
  "order": 2
}
