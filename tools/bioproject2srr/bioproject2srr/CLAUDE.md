# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This directory contains two complementary Galaxy tools for working with NCBI metadata:

### bio2srr (BioProject → SRR)
Galaxy tool wrapper for `bio2srr.py`, which retrieves SRR (Sequence Read Archive) accessions and biosample metadata from NCBI BioProjects using the Entrez E-utilities API. The tool recursively follows links to subprojects and biosamples, producing two outputs:
- `accessions.txt`: deduplicated list of SRR accessions
- `metadata.tsv`: full metadata table joining biosample attributes with SRA run information

### srr2sam (SRR → BioSample)
Galaxy tool wrapper for `srr2sam.py`, which performs the reverse operation: given a list of SRR accessions, retrieves associated BioSample metadata. Produces:
- `output.tsv`: tabular file with SRR accession as first column, followed by all BioSample metadata fields

Both tools are containerized and run via Docker in Galaxy with the image `docker.io/crashfrog/bp2srr-galaxy:latest`.

## Key Architecture

### Core Logic Flow (bio2srr.py)
1. **ID Resolution**: Convert BioProject accession → UID via esearch
2. **Recursive Link Following**: `resolve_bioproject_ids_and_links()` walks BioProject→BioProject and BioProject→BioSample links
3. **XML Parsing**: `flatten_biosample_xml()` extracts structured metadata from BioSample XML
4. **SRA Mapping**: `yield_sra_runs_from_sample()` links BioSamples to SRA runs via elink
5. **Output Generation**: Produces sorted, deduplicated TSV with custom header ordering

### E-utilities API Pattern
All NCBI API calls use `requests.get()` with:
- Rate limiting: 1s delay (no API key) or 0.1s (with `NCBI_API_KEY`)
- Batching: process IDs in groups of 200 for esummary calls
- Error handling: `raise_for_status()` on all responses

### Core Logic Flow (srr2sam.py)
1. **SRR → SRA UID**: Convert SRR accession → SRA UID via esearch
2. **SRA → BioSample Link**: Use elink to traverse from SRA to BioSample database
3. **BioSample Metadata Fetch**: Retrieve BioSample XML via esummary
4. **XML Parsing**: Reuse `flatten_biosample_xml()` from bio2srr.py
5. **Output Generation**: TSV with SRR accession as first column

### Header Ordering
Both tools use `header_sort_override()` to enforce priority fields at the left:
- **bio2srr**: `bioproject, srr_accession, biosample_accession, organism, taxid, package` then alphabetical
- **srr2sam**: `sra_run_accession, biosample_accession, organism, taxid, package` then alphabetical

## Testing

**Run unit tests:**
```bash
pytest tests.py -v
```

**Run Galaxy tool tests:**
Both tools include functional tests that are run by Galaxy's Planemo test framework. Test data is in `test-data/`.

**bio2srr.xml:**
- Input: `PRJNA681235` (BioProject accession)
- Expected outputs: `accessions.txt`, `metadata.tsv`
- Negative test: `NOTHING` (expects failure)

**srr2sam.xml:**
- Input: `srr_test_input.txt` (3 SRR accessions)
- Expected output: `biosample_metadata.tsv`

To run Galaxy tool tests locally:
```bash
planemo test --galaxy_python_version 3.11 bio2srr.xml
planemo test --galaxy_python_version 3.11 srr2sam.xml
```

## Tool Shed Publishing

This tool is published to the Galaxy ToolShed with metadata in `.shed.yml`:
- Owner: jpayne
- Category: Tools
- Name: bioproject2srr

## Dependencies

The Python script uses:
- `requests` for HTTP calls to NCBI E-utilities
- `xml.etree.ElementTree` for parsing BioSample and SRA run XML
- `itertools.batched` (or backport for Python <3.12)
- Environment variable `NCBI_API_KEY` (optional, increases rate limit)

## Important Constraints

**bio2srr:**
- Handles unreleased samples (samples without SRA runs) by including them in metadata but warning about omission from accessions list
- Duplicate SRA runs (reached via multiple projects/samples) are deduplicated in `accessions.txt` but retained in metadata
- Non-public SRA runs are skipped with a warning
- API responses contain "deliriously" embedded XML fragments that require wrapping in `<data>` root tags before parsing

**srr2sam:**
- If an SRR accession cannot be resolved to a BioSample, the row will contain "NOT_FOUND" in the biosample_accession column
- If BioSample XML parsing fails, the row will contain "PARSE_ERROR"
- Takes the first BioSample UID when an SRA record links to multiple BioSamples (rare edge case)
