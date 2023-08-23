#name: Julia Date
#description: datetime input/output
#language: julia
#input: datetime input_datetime
#output: datetime output_datetime

output_datetime = input_datetime + Dates.Day(1)
