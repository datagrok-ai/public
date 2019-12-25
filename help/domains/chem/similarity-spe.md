<!-- TITLE: Molecular Similarity Analysis using Stochastic Proximity Embedding (SPE) -->
<!-- SUBTITLE: -->

# Molecular Similarity Analysis Using Stochastic Proximity Embedding (SPE)

Molecular similarity refers to the similarity of chemical elements, molecules or chemical compounds 
with respect to either structural or functional qualities, i.e. the effect that the chemical compound 
has on reaction partners in inorganic or biological settings. Biological effects and thus also 
similarity of effects are usually quantified using the biological activity of a compound. 
In general terms, function can be related to the chemical activity of compounds (among others).

The analysis uses molecular fingerprints calculated using [Molecular Fingerprints](fingerprints.md) dialog 
and converting cross-similarities into 2D or 3D coordinates using Stochastic Proximity Embedding (SPE) 
algorithm. 

SPE is a self-organizing algorithm for producing meaningful underlying dimensions from proximity data. 
It attempts to generate low-dimensional Euclidean embeddings that best preserve the similarities 
between a set of related molecular similarities.

Similarity between pair of molecules can be estimated by different cores. 

## Available distance metrics

  * Tanimoto
  * Dice
  * Cosine
  * Sokal
  * Russel
  * Rogot/Goldberg
  * Kulczynski
  * McConnaughey
  * Asymmetric
  * Braun/Blanquet

See also:

  * [Chemical Similarity](https://en.wikipedia.org/wiki/Chemical_similarity)
  * [Molecular Fingerprints](fingerprints.md)
  * [Fingerprints - Screening and Similarity](http://www.daylight.com/dayhtml/doc/theory/theory.finger.html)
  * [Stochastic Proximity Embedding](https://pdfs.semanticscholar.org/aeb7/aa3b9655838e00de12e33e64f9f1b43bb922.pdf)
