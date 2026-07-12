# NLP changelog

## v.next

* GROK-18695: Migrated AWS Translate to AWS SDK v3 (removes unpatchable aws-sdk v2 advisory)
* GROK-18695: Dependency security updates — overrode uuid to 11.x and protobufjs to 7.6.5 (clears aws-sdk and @xenova/transformers transitive advisories)
* Docker: Cleared reported CVEs — added `apt upgrade` for base-image OS packages, upgraded pip/setuptools/wheel, and added security floors for pillow/urllib3/requests/werkzeug/filelock/idna

## 1.2.4 (2025-07-04)

Updated default params for text clustering

## 1.0.9 (2024-04-05)

Fixed multiple filters creating.

## 1.0.8 (2024-02-23)

Updated help: stemming-based search tools description is added.

## 1.0.7 (2024-02-20)

Stemming-based search (SBS) tools for text columns are added:

* SBS context & distance edit panels
* SBS approach for text embeddings computation

## 1.0.6 (2023-07-24)

Initial release of NLP package
