<!-- TITLE: Exercises -->
<!-- SUBTITLE: -->

# Exercises

This is a set of programming exercises designed make developers proficient with the
Datagrok platform. The exercises are organized as progressive steps, with tasks
of increasing complexity building on top of the previously completed steps.

During this course, we will be building support for
handling DNA nucleotide sequences. Let that not scare you, think of them as regular
strings that can only contain characters G, A, C, and T (and now you know the origins
of the "Gattaca" movie name). We will start with writing a stand-alone functions, then
with automatically recognizing nucleotide sequences in the imported data, and then 
going all the way to custom visualizations, querying relational databases, 
predictive models, integration with the external utilities, data augmentation and 
custom applications. 

## Setting up the environment

Pre-requisites: very basic JavaScript knowledge

1. Install the necessary tools
2. Get dev key for https://dev.datagrok.ai (you will work with this server)   
3. Create a default package called <Name>Sequence using datagrok-tools
4. Upload it to the server
5. Run the platform and run the package's "test" using different methods: 
    * via the Functions view
    * via the package ('Content' pane in the property panel)
    * via the console  

## Semantic types

Pre-requisites: very basic JavaScript knowledge

1. Create a `complement` function that takes a nucleotide string and returns it complement.
   Essentially, change each character to the complementary one: A<=>T, G<=>C. 
   Run it and check whether everything works fine. 
2. Now, let's specify that this function is meant to accept not any string, but nucleotides only,
   and returns a nucleotide string as well. In order to do that, let's annotate both input and output parameters with the `dna_nucleotide`
   semantic type. At this point, `dna_nucleotide` string does not have any meaning, but we'll connect
   the dots later.
3. Define a `detectNucleotides` semantic type detector function as part of the special
   `detectors.js` file. It should check whether a column is a string column, and whether
   each string represents a nucleotide (hint: for best performance, 
   don't iterate over all column values, instead iterate on categories). When everything is done
   correctly, the `detectors.js` file will get loaded by the platform automatically, and the
   `detectNucleotides` function will be executed against every column in a newly added table.
4. Test your implementation by opening the following CSV file (or go Open | Text, and paste it there):
   ```
   sequence, id
   GATTACA, 1997
   ATTCGGA, 1984
   TTTAGGC, 2021 
   ```
   Hover over the `sequence` column header after the data is imported - if everything is done correctly,
   you will see 'quality: dna_nucleotide' in the bottom of the tooltip. Alternatively, you can 
   find this information if you click on the column and expand the "Details" pane on the right.
5. Now tag the previously created `complement` function with a `panel` tag. This will instruct the
   platform to use the `complement` function for providing additional information for string values
   of the `dna_nucleotide` semantic type. To test it, simply open our test file, click on any cell
   in the `sequence` column, and find the `complement` property panel on the right.