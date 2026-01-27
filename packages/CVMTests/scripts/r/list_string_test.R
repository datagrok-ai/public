#name: R List String Test
#description: Returns the last element of a string list
#language: r
#input: list<string> string_list
#output: string result

result <- string_list[length(string_list)]
