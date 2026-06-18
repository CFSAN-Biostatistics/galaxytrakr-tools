"Retrieve BioSample metadata from SRR accessions via Eutils"

import requests
import sys
import csv
import os

try:
    from itertools import batched
except ImportError:
    from itertools import islice
    def batched(iterable, n):
        "Batch data into tuples of length n. The last batch may be shorter."
        if n < 1:
            raise ValueError('n must be at least one')
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch

from functools import cmp_to_key
from time import sleep
from xml.etree import ElementTree as xml

esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
elink = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

import logging
logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("srr2sam")

extra_params = {}

api_key = os.environ.get("NCBI_API_KEY")

if api_key:
    logger.info(f"Using NCBI API key {api_key[:4]}{'*' * (len(api_key) - 8)}{api_key[-4:]}")
    extra_params["api_key"] = api_key

def log(msg):
    if api_key:
        logger.info(msg.replace(api_key, f"{api_key[:4]}{'*' * (len(api_key) - 8)}{api_key[-4:]}"))
    else:
        logger.info(msg)

def get_tag(root, tag):
    val = root.find(tag)
    if val is not None:
        return val.text
    log(f"No result for {tag}")

def header_sort_override(a, b):
    if a == b:
        return 0
    try:
        for name in ["sra_run_accession", "biosample_accession", "organism", "taxid", "package"]:
            if a == name:
                return -1
            if b == name:
                return 1
    except:
        pass
    if a < b:
        return -1
    else:
        return 1

hso = cmp_to_key(header_sort_override)

def flatten_biosample_xml(biosampxml):
    "Parse BioSample XML and return dict of attributes"
    root = xml.fromstring(biosampxml)
    accession = get_tag(root, r'.//Id[@db="BioSample"]')
    organism = get_tag(root, r".//OrganismName")
    tax_id = root.find(r".//Organism").attrib.get("taxonomy_id")
    package = get_tag(root, r".//Package")
    sampledict = dict(
        biosample_accession=accession,
        organism=organism,
        taxid=tax_id,
        package=package
    )
    for attribute in root.findall("Attributes/Attribute"):
        sampledict[attribute.attrib.get("harmonized_name", attribute.attrib['attribute_name'])] = attribute.text

    return sampledict

def srr_to_biosample_uid(srr_accession):
    "Convert SRR accession to BioSample UID via Entrez"
    sleep(1 if not api_key else 0.1)

    # First convert SRR accession to SRA UID
    response = requests.get(esearch, params=dict(db="sra", term=srr_accession, format="json", **extra_params))
    response.raise_for_status()
    reply = response.json()

    try:
        sra_uid = reply["esearchresult"]["idlist"][0]
        log(f"Found SRA UID {sra_uid} for '{srr_accession}'")
    except IndexError:
        logger.warning(f"No SRA UID found for '{srr_accession}'")
        return None

    sleep(1 if not api_key else 0.1)

    # Link from SRA to BioSample
    response = requests.get(elink, params=dict(id=sra_uid, dbfrom="sra", db="biosample", format="json", **extra_params))
    response.raise_for_status()
    reply = response.json()

    biosample_uids = reply.get("linksets", [{}])[0].get("linksetdbs", [{}])[0].get("links", [])

    if not biosample_uids:
        logger.warning(f"No BioSample link found for SRR {srr_accession} (UID {sra_uid})")
        return None

    return biosample_uids[0]  # Return first BioSample UID

def get_biosample_metadata(biosample_uid):
    "Fetch BioSample metadata XML via esummary"
    sleep(1 if not api_key else 0.1)

    response = requests.get(esummary, params=dict(id=biosample_uid, db="biosample", format="json", **extra_params))
    response.raise_for_status()
    reply = response.json()

    biosample_xml = reply.get("result", {}).get(biosample_uid, {}).get("sampledata")

    if not biosample_xml:
        logger.warning(f"No metadata found for BioSample UID {biosample_uid}")
        return None

    return biosample_xml

def main(input_file):
    rows = []

    # Read SRR accessions from input file
    with open(input_file, 'r') as f:
        srr_accessions = [line.strip() for line in f if line.strip()]

    log(f"Processing {len(srr_accessions)} SRR accessions")

    for i, srr in enumerate(srr_accessions, 1):
        log(f"Processing {srr} ({i}/{len(srr_accessions)})")

        # Get BioSample UID
        biosample_uid = srr_to_biosample_uid(srr)
        if not biosample_uid:
            # Create minimal row for failed lookups
            rows.append({"sra_run_accession": srr, "biosample_accession": "NOT_FOUND"})
            continue

        # Get BioSample metadata
        biosample_xml = get_biosample_metadata(biosample_uid)
        if not biosample_xml:
            rows.append({"sra_run_accession": srr, "biosample_accession": "NOT_FOUND"})
            continue

        # Parse metadata
        try:
            metadata = flatten_biosample_xml(biosample_xml)
            metadata["sra_run_accession"] = srr
            rows.append(metadata)
        except Exception as e:
            logger.error(f"Error parsing BioSample metadata for {srr}: {e}")
            rows.append({"sra_run_accession": srr, "biosample_accession": "PARSE_ERROR"})

    log(f"Writing {len(rows)} rows to output.tsv")

    # Collect all headers
    header = set()
    for row in rows:
        for key in row.keys():
            header.add(key)

    header = sorted(list(header), key=hso)

    # Sort by SRR accession
    rows.sort(key=lambda x: x.get("sra_run_accession", ""))

    # Write output
    with open("output.tsv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", dialect="excel")
        writer.writeheader()
        writer.writerows(rows)

    log("Complete")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        logger.error("Usage: python srr2sam.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.exists(input_file):
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    try:
        main(input_file)
    except requests.HTTPError as e:
        logger.error(f"HTTP error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)
