#name: Octave Params Test
#language: octave
#tags: test, selenium
#input: bool BOOL_INPUT = true
#input: string STRING_INPUT = input_string
#input: double DOUBLE_INPUT = 3.1415926
#input: int INT_INPUT = 123456
#input: dataframe table [Data table]
#input: column COLUMN_INPUT
#input: column_list COLUMN_LIST_INPUT
#output: bool BOOL_OUTPUT
#output: string STRING_OUTPUT
#output: double DOUBLE_OUTPUT
#output: int INT_OUTPUT
#output: dataframe DF_OUTPUT
BOOL_OUTPUT = ~BOOL_INPUT
STRING_OUTPUT = [STRING_INPUT " -" COLUMN_LIST_INPUT{1}]
DOUBLE_OUTPUT = 2 * DOUBLE_INPUT
INT_OUTPUT = rows(table) * columns(table) + INT_INPUT
DF_OUTPUT = {1, "first"; 3, "second"}