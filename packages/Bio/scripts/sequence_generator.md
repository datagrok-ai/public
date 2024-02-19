# Sequence generator

The utility generates clusters of macromolecule sequences to test SAR functionality. 
Each cluster contains a randomly generated sequence motif.
Each sequence has activity - a Gauss-distributed random value. 
The utility can simulate activity cliffs - random changes in the conservative motif letters,
leading to the significant change in the activity.
Utility can simulate multiple experimental assays measuring activity, with different scales and noise levels.

### Run options
The utility can work in two modes:
* Standalone command-line tool. Run the utility with `--help` key to get detailed help
* Datagrok script. In this mode, Datagrok automatically generates utility UI.

## Utility algorithm

### Motif generation
* Specify the sequence alphabet: DNA/RNA/Peptides/HELM file
* Generate motif template
  * Iterate through all positions in the motif
  * With the probability `prob_any` it will be a non-conservative letter.
    In this case, place the `?` code in this position.
    The non-conservative letters can't reside on motif ends.
  * Otherwise, generate a conservative letter
    * Generate a random integer in range `[1-max_variants_cluster]` - the number of possible letters in this position.
    * Randomly choose corresponding of letters from the alphabet.  
* Generate cluster sequences by template
  * Iterate through all motif template positions
  * Take the letter from the list of allowed letters for the position
* With the probability `cliff_probability` generate sequence cliff:
  * Choose a conservative letters
  * Mutate it to a random letter not present in the allowed template letters for this position.

## Generation activities 

The activity generation simulates experimental data from some real assay with noise and bias.
* Generate random mean value for Gauss distribution for "ideal" activity generation. 
    The `activity_range` parameter defines the range of ideal activities
* Generate "ideal" activities from gauss distribution with sigma=1
* Make activity cliffs for the of sequences:
  * Calculate difference between generated activities
  * Calculate `cliff_size` - random Gauss-distributed value   
  * "Push apart" these values to ensure that the difference between them is 
* Calculate "assay activities" for each assay
  * Generate "noise" - Gauss-distributed value with a randomly chosen center    
  * Calculate `real_activity = ideal_activity + noise*noise_level`
    The `noise_level` is the ration of signal/noise for the assay
  * Rescale the activity to fit in the requested assay scale 
