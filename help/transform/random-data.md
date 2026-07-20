---
title: "Random data"
description: Add a numerical column filled with random data from a chosen probability distribution and seed.
keywords:
  - generate random data
  - random column
  - synthetic dataset
  - normal distribution sample
  - random number generation
  - test data generator
---

Adds a numerical column with random data with the specified distribution with initial seed. Parameters of the
distribution can be edited as well.

This can be useful for hypothesis testing, modeling, as well as for generating synthetic datasets.

Supported distributions and parameters:

| Distribution | Parameters |
|--------------|------------|
| normal       | mean, sd   |
| log-normal   | mean, sd   |
| binomial     | size, prob |
| poisson      | lambda     |
| uniform      | min, max   |

See also:

* [Randomness](https://en.wikipedia.org/wiki/Randomness)
* [List of probability distributions](https://en.wikipedia.org/wiki/List_of_probability_distributions)
* JS API: [Random data](https://public.datagrok.ai/js/samples/domains/data-science/random-data)
