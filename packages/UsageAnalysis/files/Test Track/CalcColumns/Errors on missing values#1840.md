Data: [ForDGWithSigns.csv](https://github.com/datagrok-ai/public/files/11324437/ForDGWithSigns.csv)

1. Add new columns with the following formulas:

    * Qnum((100-${withSigns2})/100, if(qualifier(${withSigns2})==">", "<", "="))
    * Qnum((100-${withSigns2})/100, if(isnotempty("${withSigns2}"), if(qualifier(${withSigns2})==">", "<", "="), null))

2. Check that no errors were genrated by The Qnum() function, especially on missing values

---

{
  "order": 2
}
