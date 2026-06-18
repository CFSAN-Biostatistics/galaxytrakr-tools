"Grab SRR numbers from Bioprojects and sub-bioprojects via Eutils"

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
        # batched('ABCDEFG', 3) --> ABC DEF G
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

logger = logging.getLogger("bio2srr")

extra_params = {}

api_key = os.environ.get("NCBI_API_KEY")

if api_key:
    logger.info(f"Using NCBI API key {api_key[:4]}{'*' * (len(api_key) - 8)}{api_key[-4:]}")
    extra_params["api_key"] = api_key

def log(msg):
    if api_key:
        logger.info(msg.replace(api_key, f"{api_key[:4]}{'*' * (len(api_key) - 8)}{api_key[-4:]}")) # fix logging later
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
        for name in ["bioproject", "srr_accession", "biosample_accession", "organism", "taxid", "package",]:
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

def resolve_bioproject_ids_and_links(bioproject_id_list):
    "Recursively follow bioproject and biosample links, yield biosample UID's and biosample XML"
    for i, (bioproject, bioproject_id) in enumerate(bioproject_id_list):
        log(f"Processing {bioproject} ({bioproject_id}) {i+1}/{len(bioproject_id_list)}")
        #get bioproject to bioproject links
        response = requests.get(elink, params=dict(db="bioproject", dbfrom="bioproject", id=bioproject_id, format="json", **extra_params))
        response.raise_for_status()
        reply = response.json()
        linksets = reply.get("linksets", [{}])[0].get("linksetdbs", [0,0,{}])
        if len(linksets) >= 3:
            for id in linksets[2].get("links", []): #third index is the up to down links
                response = requests.get(esummary, params=dict(id=id, db="bioproject", format="json"))
                response.raise_for_status()
                replyy = response.json()
                biop = replyy["result"][id]["project_acc"]
                if id not in bioproject_id_list:
                    bioproject_id_list.append((biop, id)) # recurse over bioproject links
        # get bioproject to biosample links
        response = requests.get(elink, params=dict(db="biosample", dbfrom="bioproject", id=bioproject_id, format="json", **extra_params))
        response.raise_for_status()
        reply = response.json()
        links = reply.get("linksets", [{}])[0].get("linksetdbs", [{}])[0].get("links", [])
        log(f"Found {len(links)} biosample links for {bioproject} ({bioproject_id})")
        for ids in batched(links, 200):
            response = requests.get(esummary, params=dict(id=",".join(ids), db="biosample", format="json"))
            response.raise_for_status()
            replyy = response.json()
            for field, value in replyy.get("result", {}).items():
                if "uids" not in field:
                    yield bioproject, field, value["sampledata"] # this is XML, deleriously
                    sleep(1 if not api_key else 0.1)


biosample_example = """
<BioSample access="public" publication_date="2020-12-21T00:00:00.000" last_update="2022-06-23T17:45:35.674" submission_date="2020-12-21T15:08:05.690" id="17131268" accession="SAMN17131268">
    <Ids>     
        <Id db="BioSample" is_primary="1">SAMN17131268</Id>
        <Id db_label="Sample name">CJP19-D996</Id>
    </Ids>
    <Description>     
        <Title>Pathogen: environmental/food/other sample from Campylobacter jejuni</Title>     
        <Organism taxonomy_id="197" taxonomy_name="Campylobacter jejuni">       
            <OrganismName>Campylobacter jejuni</OrganismName>
        </Organism>   
    </Description>   
    <Owner>     
        <Name url="http://www.fda.gov/Food/FoodScienceResearch/WholeGenomeSequencingProgramWGS/default.htm" abbreviation="CFSAN">FDA Center for Food Safety and Applied Nutrition</Name> 
    </Owner> 
    <Models>   
        <Model>Pathogen.env</Model>  
    </Models>  
    <Package display_name="Pathogen: environmental/food/other; version 1.0">Pathogen.env.1.0</Package> 
    <Attributes>  
        <Attribute attribute_name="strain" harmonized_name="strain" display_name="strain">CJP19-D996</Attribute>  
        <Attribute attribute_name="collection_date" harmonized_name="collection_date" display_name="collection date">missing</Attribute>     
        <Attribute attribute_name="geo_loc_name" harmonized_name="geo_loc_name" display_name="geographic location">missing</Attribute>     
        <Attribute attribute_name="collected_by" harmonized_name="collected_by" display_name="collected by">CDC</Attribute>     
        <Attribute attribute_name="lat_lon" harmonized_name="lat_lon" display_name="latitude and longitude">missing</Attribute>     
        <Attribute attribute_name="isolation_source" harmonized_name="isolation_source" display_name="isolation source">missing</Attribute>     
        <Attribute attribute_name="isolate" harmonized_name="isolate" display_name="isolate">CFSAN091032</Attribute>     
        <Attribute attribute_name="project name" harmonized_name="project_name" display_name="project name">GenomeTrakr</Attribute>     
        <Attribute attribute_name="sequenced by" harmonized_name="sequenced_by" display_name="sequenced by">FDA Center for Food Safety and Applied Nutrition</Attribute>   
    </Attributes>   
    <Links>     
        <Link type="entrez" target="bioproject" label="PRJNA681235">681235</Link>   
    </Links>
    <Status status="live" when="2020-12-21T15:08:05.693"/> 
</BioSample>

"""

def flatten_biosample_xml(biosampxml):
    root = xml.fromstring(biosampxml)
    accession = get_tag(root, r'.//Id[@db="BioSample"]')
    # sample_name = get_tag(root, r'.//Id[@db_label="Sample name"]')
    organism = get_tag(root, r".//OrganismName")
    tax_id = root.find(r".//Organism").attrib.get("taxonomy_id")
    package = get_tag(root, r".//Package")
    sampledict = dict(
        biosample_accession=accession,
        # sample_name=sample_name,
        organism = organism,
        taxid = tax_id,
        package = package
    )
    for attribute in root.findall("Attributes/Attribute"):
        sampledict[attribute.attrib.get("harmonized_name", attribute.attrib['attribute_name'])] = attribute.text

    return sampledict


def yield_sra_runs_from_sample(biosample):
    sleep(1 if not api_key else 0.1)
    response = requests.get(elink, params=dict(id=biosample, dbfrom="biosample", db="sra", format="json", **extra_params))
    response.raise_for_status()
    reply = response.json()
    for ids in batched(reply.get("linksets", [{}])[0].get("linksetdbs", [{}])[0].get("links", []), 200):
        sleep(1 if not api_key else 0.1)
        response = requests.get(esummary, params=dict(id=','.join(ids), db="sra", format="json", **extra_params))
        response.raise_for_status()
        replyy = response.json()
        for field, value in replyy.get("result", {}).items():
            if "uids" not in field:
                yield field, value.get("runs")


runs_example = """
<Run acc="SRR13167188" total_spots="827691" total_bases="385043067" load_done="true" is_public="true" cluster_name="public" static_data_available="true"/>  
<Run acc="SRR13167189" total_spots="827691" total_bases="385043067" load_done="true" is_public="true" cluster_name="public" static_data_available="true"/>   
"""

def flatten_runs(runxml):
    root = xml.fromstring(f"<data>{runxml}</data>") # gotta fix their garbage embedded XML since it isn't singly-rooted
    for run in root.findall(".//Run"):
        if run.attrib["is_public"] == "false":
            logger.warning(f"Skipping non-public run {run.attrib['acc']}")
        yield dict(
            sra_run_accession = run.attrib["acc"],
            total_spots = run.attrib["total_spots"],
            total_bases = run.attrib["total_bases"],
        )



def main(starting_bioproject):
    rows = []
    response = requests.get(esearch, params=dict(db="bioproject", term=starting_bioproject, field="PRJA", format="json"))
    response.raise_for_status()
    reply = response.json()
    try:
        bioproject_id = reply["esearchresult"]["idlist"][0]
        log(f"Found UID {bioproject_id} for '{starting_bioproject}'")
    except IndexError:
        logger.error(f"No results found for '{starting_bioproject}'. Error was \"{reply['esearchresult']['warninglist']['outputmessages']}\"")
        sys.exit(1)
    sleep(1 if not api_key else 0.1)
    for bioproject, biosample, biosample_xml in resolve_bioproject_ids_and_links([(starting_bioproject, bioproject_id)]):
        try:
            sampledict = flatten_biosample_xml(biosample_xml)
        except KeyError:
            log(biosample_xml)
            raise
        sampledict["bioproject"] = bioproject
        noruns = True
        for sra, runs in yield_sra_runs_from_sample(biosample):
            for run in flatten_runs(runs.strip()):
                noruns = False
                run.update(sampledict)
                rows.append(run)
        if noruns:
            rows.append(sampledict)

    log(f"Writing {len(rows)} rows to metadata.tsv")

    header = set()
    for row in rows:
        for key in row.keys():
            header.add(key)

    header = sorted(list(header), key=hso)
    # logger.info(f"Header: {header}")

    rows.sort(key=lambda x: x["biosample_accession"])

    with open("metadata.tsv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", dialect="excel")
        writer.writeheader()
        writer.writerows(rows)

    # check for duplicate runs and unreleased samples

    accessions = [row.get("sra_run_accession") for row in rows if row.get("sra_run_accession")]

    raw_length = len(accessions)

    accessions = sorted(list(set(accessions)))

    if raw_length < len(rows):
        logger.warning(f"Bioproject {starting_bioproject} contains unreleased samples. {len(rows) - raw_length} samples will not be included in accessions.txt")

    if len(accessions) < raw_length:
        logger.warning(f"Some SRA runs may have been reached through multiple projects or samples. accessions.txt will be deduplicated but the metadata table is not")

    log(f"Writing {len(accessions)} unique accessions to accessions.txt")

    with open("accessions.txt", "w") as f:
        for accession in accessions:
            f.write(accession + "\n")


if __name__ == "__main__":
    b = sys.argv[1].strip()
    log(f"Starting with {b}")
    try:
        main(b)
    except requests.HTTPError as e:
        logger.error(e)
        sys.exit(1)

                



