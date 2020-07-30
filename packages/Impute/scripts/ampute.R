#name: ampute
#description: simulate missing values
#help-url: https://www.rdocumentation.org/packages/mice/versions/3.9.0/topics/ampute
#language: r
#tags: demo, viewers
#input: dataframe data [Input data table]
#input: column_list columns [list of all columns of interest]
#input: double prob = 0.2 [NA proportion]
#input: string mechanism {choices: ["MCAR", "MAR","MNAR"]} [missingness mechanism]
#output: dataframe NA_sim_df [df with simulated missing values]
require(mice)

data <- data[,columns]
data <- data[,colSums(is.na(data))<nrow(data)]

results <- ampute(data,
                  prop = prob,
                  patterns = NULL,
                  freq = NULL,
                  mech = mechanism,
                  weights = NULL,
                  std = TRUE,
                  cont = TRUE,
                  type = NULL,
                  odds = NULL,
                  bycases = TRUE,
                  run = TRUE)

NA_sim_df  <- results$amp