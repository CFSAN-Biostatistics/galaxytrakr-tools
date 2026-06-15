#!/usr/bin/env python3

#
# Script to retrieve shed.yml metadata fields from a
# galaxy tool repository (opentrakr.org). These can
# then be used with planemo to push tools to new repos
# Errol Strain (estrain@gmail.com)
#

import os
import re
import textwrap
import requests
import yaml
from bs4 import BeautifulSoup

owner = "estrain"
tools_dir = "tools"
base_html = "https://opentrakr.org/repository/view_repository?id="

requests.packages.urllib3.disable_warnings()

def fetch_repo_info(tool_name):
    api = f"https://opentrakr.org/api/repositories?owner={owner}&name={tool_name}"
    r = requests.get(api, verify=True)
    if not r.ok:
        return None
    data = r.json()
    if not data:
        return None
    return data[0]

def clean_text(t):
    if not t:
        return ""
    t = t.replace("\xa0", " ")     # non-breaking space
    t = re.sub(r"\s+", " ", t)
    return t.strip()

def parse_html_description(repo_id):
    html_url = f"{base_html}{repo_id}"
    r_html = requests.get(html_url, verify=True)
    if not r_html.ok:
        return "", "", None

    soup = BeautifulSoup(r_html.text, "html.parser")

    short_desc = ""
    long_desc = ""
    homepage = None

    # Synopsis
    syn_div = soup.find("b", string=lambda s: s and "Synopsis" in s)
    if syn_div:
        if syn_div.next_sibling and isinstance(syn_div.next_sibling, str):
            short_desc = clean_text(syn_div.next_sibling)

    # Detailed description
    label = soup.find("label", string=lambda s: s and "Detailed description" in s)
    if label:
        desc_table = label.find_next("table", id="description_table")
        if desc_table:
            td = desc_table.find("td")
            if td:
                long_desc = clean_text(td.get_text(" ", strip=True))

    # Homepage
    for a in soup.find_all("a", href=True):
        if "github" in a["href"].lower():
            homepage = a["href"]
            break

    return short_desc, long_desc, homepage

class LiteralStr(str): pass
def literal_str_representer(dumper, value):
    return dumper.represent_scalar('tag:yaml.org,2002:str', value, style='|')
yaml.add_representer(LiteralStr, literal_str_representer)

for tool in sorted(os.listdir(tools_dir)):
    tool_path = os.path.join(tools_dir, tool)
    if not os.path.isdir(tool_path):
        continue

    shed_file = os.path.join(tool_path, ".shed.yml")
    if os.path.exists(shed_file):
        print(f"Skipping existing {shed_file}")
        continue

    repo_info = fetch_repo_info(tool)
    if not repo_info:
        print(f"No repository found for {tool}")
        continue

    repo_id = repo_info.get("id")
    short_desc, long_desc, homepage = parse_html_description(repo_id)

    # Wrap long description for readability
    wrapped = "\n".join(textwrap.wrap(long_desc, width=80)) if long_desc else ""

    data = {
        "name": tool,
        "owner": owner,
        "description": short_desc or repo_info.get("description", ""),
        "long_description": LiteralStr(wrapped),
        "homepage_url": homepage or repo_info.get("homepage_url", ""),
        "type": repo_info.get("type", "unrestricted"),
        "categories": [],
    }

    with open(shed_file, "w") as fh:
        yaml.dump(data, fh, sort_keys=False)
    print(f"Created {shed_file}")
