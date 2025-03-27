---
title: "Performance tips"
---

Datagrok was designed to be as efficient as possible. While many of the features work very fast any way you use them,
others either need a deeper understanding of the internal structures, or need to be fine-tuned to unlock maximum
performance. This document contains a set of guiding principles and approaches targeted for application programmers.

as an application programmer

## DataFrame

DataFrame is a [highly optimized columnar data engine](../under-the-hood/scaling.md#in-memory-database) that is used for data ingestion and
transfer, data transformations, visualization, and analyses. To understand the

* For maximum performance, use Column.getRawData() method

## Databases

* Caching

## Static datasets

## UI

Datagrok's UI should be friendly and snappy.
