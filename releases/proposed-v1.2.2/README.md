# Proposed FTL methods v1.2.2: owner review kit

Prepared 2026-10-02. **This is a proposal, not a published Zenodo release.** The proposed record type is publication/report. The latest public FTL release remains [v1.2.1](https://doi.org/10.5281/zenodo.20218030), in [FTL concept 17726159](https://doi.org/10.5281/zenodo.17726159).

Download the [methods ZIP](FTL_METHODS_v1_2_2_20261002.zip) and [proposed Zenodo metadata](ZENODO_LEGACY_METADATA_PROPOSED.json). The JSON preserves the prepared legacy metadata exactly; its fields are inside the top-level `metadata` object. It contains no authentication information.

The report collects reply-delay and unequal-speed calculations, an explicitly assumed preferred-time free scalar exercise, and their reproducible internal checks. These are conditional mathematical and software results. This is not empirical v1.3.0, a demonstrated FTL channel, a novelty claim, or external peer review.

## Package and source checks

The ZIP is 169,189 bytes. Its SHA-256 is:

```text
1d80226cb9e8c30e43604ee650b076d04f6e9faebac492819df674228a0560f7
```

ZIP CRC and source-hash checks passed. The package contains 40 original public source files from [reviewed main commit 7ad6a32](https://github.com/maldonado-research/faster-than-light/tree/7ad6a32a2f413b8f13c329ae685478b19ea42dd2), plus five release cover/provenance files. `SOURCE_MANIFEST.json` inside the ZIP binds those 40 source files. The originals are byte-identical to that commit.

The historical `CITATION.cff`, September research/documentation snapshots, and dated source records are preserved. `RELEASE_CITATION.cff` is a separate provisional citation for the proposed edition and asserts no new version DOI. Included registrations, results and internal audits describe the scope of the computational checks; they are not physical observations.

The package and this review folder contain no private raw archives, transcripts, downloaded article contents, original downloaded toolkits, observational posterior bundles, credentials, or authentication/error logs. Documentation and generated data retain CC BY 4.0; original code retains MIT under the packaged `LICENSE`.

## Finish the owner draft

The [owner draft 23111786](https://zenodo.org/uploads/23111786) is incomplete: its metadata is empty and it contains two inherited toolkit copies. It is not ready to publish.

1. Download and inspect the ZIP and proposed metadata. Confirm the proposed v1.2.2 version, report category, author, description, licence and actual publication date.
2. Open the draft while signed in to the owning Zenodo account. Fill its fields from the prepared metadata; the JSON is a field reference, not an authentication file.
3. In **this draft only**, upload this methods ZIP. Verify its filename, 169,189-byte size and checksum against the [file reference](ZENODO_FORM_FIELDS.md#verify-the-uploaded-package), then remove only the two inherited toolkit references. Save, reopen and verify the fields and single new file.
4. Review the complete draft, its files, claim boundaries, citation and relationship to the existing FTL concept.
5. After reviewing the complete metadata and single verified package, you can choose Publish when satisfied. Until the owner publishes and the public record is verified, v1.2.2 remains proposed.

Replacing files in this unpublished draft leaves the old public v1.2.1 record intact. No deletion or change to that public record is requested here.

## Integration repair and current API readback

The owner's later signed-in browser check supersedes the earlier enabled-integration report: GitHub release histories were empty for all eight public repositories, and automatic GitHub archiving was switched **OFF** for the seven research repositories and shared website. Reload confirmed all eight stayed off. This is owner-reported browser evidence; the research environment did not independently inspect those switches. Existing repositories, tokens, published records and valid drafts are preserved.

Keep automatic GitHub archiving off. Use the existing token/API workflow and the established FTL DOI family, reusing draft23111786. Do not re-enable archiving, create a test GitHub release, or create a second standalone record for this same work. Intentional versions in this family and genuinely separate scientific supplements are distinct from duplicate publications; similar titles alone do not justify deleting a draft or record.

Fresh authenticated reads in the research environment succeeded and confirmed this same owner-scoped draft in concept17726159. Its release metadata is still empty and its files are still the two inherited kits; the methods ZIP is absent. A reserved draft DOI or HTTP200 is not a published update. The earlier JSON, binary and multipart persistence failures remain unresolved; current evidence does not show missing credentials. Official documentation and offline nonempty-body/framing checks found no justified unchanged API retry.

Use the [readable form fields](ZENODO_FORM_FIELDS.md) with the existing signed-in draft. The frozen ZIP and JSON are unchanged, and this proposal does not include subsequent coupled-scalar calculations. After reviewing and verifying a complete draft, the owner may choose Publish through the process above.

Official references: [enable a repository](https://help.zenodo.org/docs/github/enable-repository/), [deposit API](https://developers.zenodo.org/) and [manage versions](https://help.zenodo.org/docs/deposit/manage-versions/).

The [main repository overview](../../README.md) remains the starting point for the research.
