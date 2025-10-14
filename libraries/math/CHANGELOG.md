# math changelog

## 1.2.6 (2024-12-17)

* Vector distance functions support on webGPU

## 1.2.5 (2024-12-17)

* Generations updated

## 1.2.4 (2024-12-13)

* Remove deprecated adapter info function

## 1.2.3 (2024-10-17)

Fix compatibility with newer version of wgpu types

## 1.2.2 (2024-10-16)

* Increase complexity of bit array functions

## 1.2.1 (2024-09-25)

* Expect navigator.gpu to be undefined in some browsers

## 1.1.13 (2024-09-06)

Bump Api version

## 1.1.12 (2024-08-07)

Correct NW calculation on webGPU

## 1.1.11 (2024-08-06)

Modify MCL and NW calculations on webGPU

## 1.1.10 (2024-06-07)

MMP generations on webGPU

## 1.1.9 (2024-05-27)

Add inflation factor to MCL

## 1.1.8 (2024-05-16)

Fix inconsistent KNN size in webGPU.

## 1.1.7 (2024-05-09)

* Fix webGPU description nullish value.

## 1.1.6 (2024-04-25)

Matched molecular pairs analysis paralleled via gpu.

## 1.1.5 (2024-04-25)

Use npm for gpu rules storage.

## 1.1.4 (2024-04-25)

Add webGPU MCL implementation.

## 1.1.3 (2024-04-23)

Introduce rules for GPU usage

## 1.1.2 (2024-04-15)

Fix webGPU numeric distance with 0 range

## 1.1.1 (2024-04-14)

Add checks for webGPU distances.

## 1.1.0 (2024-04-14)

### Features

* Add webGPU Sparse matrix calculation.
* Add webGPU UMAP implementation.
* Restructure API to be more modular.

## 1.0.8 (2024-04-05)

Add webGPU KNN calculation.

## 1.0.2 (2023-07-21)

A new library intended to be used for high efficiency (mostly wasm) calculations.

*Dependency: datagarok-api >= 1.15.0*

### Features

Added wasm algorithm for hierarchical clustering, which is used by the [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram) package.
