#name: R Map
#description: map input/output
#language: r
#input: map input_map
#input: string unique_key
#output: map output_map

input_map[[unique_key]] <- 'Datagrok'
output_map <- input_map
