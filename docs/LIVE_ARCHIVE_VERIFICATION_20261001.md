# Live verification of the FTL Zenodo family

Verified 2026-10-01 Pacific (2026-10-02 UTC). This later check supplements the September 30 edition's historical statement that live retrieval was unavailable. It preserves that edition and its source snapshot.

TLS-verified official API responses identify the latest published record in concept family **17726159** as **v1.2.1**, [DOI 10.5281/zenodo.20218030](https://doi.org/10.5281/zenodo.20218030), published 2026-05-15, titled *Einstein Audit for Faster-Than-Light Hidden Sectors (FTL-only): A2 Strict Protocol Freeze and Provider-Return Gate*.

The concept record and latest-version endpoint resolve to the same record. Both the versions endpoint and the all-versions concept query returned a total of four records, four results and no next-page link:

| Version | Record | Publication date |
|---|---|---|
| v1.0.0 | [17726160](https://zenodo.org/records/17726160) | 2025-11-26 |
| v1.1.0 | [17926499](https://zenodo.org/records/17926499) | 2025-12-13 |
| v1.2.0 | [18499411](https://zenodo.org/records/18499411) | 2026-02-05 |
| v1.2.1 | [20218030](https://zenodo.org/records/20218030) | 2026-05-15 |

This establishes the observed published concept family at verification time. It does not establish the absence of unrelated, private or unpublished records.

Official query endpoints:

- [Record metadata](https://zenodo.org/api/records/20218030)
- [Concept metadata](https://zenodo.org/api/records/17726159)
- [Latest version](https://zenodo.org/api/records/20218030/versions/latest)
- [Version list](https://zenodo.org/api/records/20218030/versions)
- [All-version concept query](https://zenodo.org/api/records?q=conceptrecid%3A17726159&all_versions=true)

## Recovered public artifacts

The record contains two public ZIPs. Downloaded sizes, declared MD5 checksums and ZIP CRCs passed. Computed SHA-256 values are recorded below. The v6.1 archive contains a byte-identical v6.0 payload; all 12 internal manifest entries passed their length and SHA-256 checks.

| Artifact | Bytes | Declared MD5 |
|---|---:|---|
| FTL_A2_STRICT_PROTOCOL_FREEZE_KIT_v6_0.zip | 1,502,150 | d9fb0abc2ae08c06bf2091171273191f |
| FTL_A2_ZENODO_PROTOCOL_RELEASE_DECISION_KIT_v6_1.zip | 1,463,506 | 45b332403cfd5b9688e5bccec27f5a24 |

Computed SHA-256:

```text
v6.0: e865df05a1aa01322020807dd31b1d835edecff26df1960f72734ca3fa3f23c3
v6.1: d8f5db7afe74db2d2419f5d503416c897fa54038c19add9a81a102be7dbba40b
```

No v6.3 kit was identified in this public family's file metadata. Public v6.1 expressly excludes strict real-data A2 claims. Its included fixture metadata labels the input as synthetic smoke-test material. File names containing “observational” do not change that provenance. The contents of all older nested archives and PDFs were not comprehensively reviewed.

## Current synthetic smoke run

The unchanged recovered v6.0 runner was executed outside the Git checkout using its included synthetic fixture. Static review checked script calls and ZIP paths before execution. The run used 150 radicand draws, 12 PCI bootstrap draws, 800 samples per class, three folds, base seed 60, W_cap=10 and a 180-second bound.

Three configurations completed: primary sigmoid/within-only, sensitivity isotonic/within-only, and sensitivity sigmoid/event-and-within. Each child exited zero and reported provider pass and COMPLETE_STABLE_5PCT. Current outputs were independently inspected rather than relying on the wrapper's exit status or archived smoke reports.

These are **synthetic software checks**. The historical classifier bias and other subsequent correction reports still apply; stability does not certify estimator accuracy or validate missing v6.3 results. This run does not establish an observational constraint, FTL detection or peer review.

The wrapper itself does not propagate child failures, and its synthetic promotion label is chosen even when a child fails. Future reproduction must inspect each child's exit/status and fresh numerical outputs. A zero wrapper exit alone is insufficient. Full parameters, environment versions and independently checked output counts are recorded in [the current smoke validation](../research/public-kit-smoke/VALIDATION.json).

No Zenodo record was edited, no new version or DOI was created, and the intended strict real-data v1.3.0 milestone remains unreached.
