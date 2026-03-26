#!/usr/bin/env python
# Python 3.6 compatible

import argparse
import json
import os
import shutil
import tarfile
import tempfile
from datetime import datetime

import wget

DB_paths = {
    "mgnify_v5_lsu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_lsu-20200130.tar.gz",
    "mgnify_v5_ssu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_ssu-20200130.tar.gz",
    "mgnify_v5_its_unite": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/UNITE-20200214.tar.gz",
    "mgnify_v5_its_itsonedb": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/ITSoneDB-20200214.tar.gz",
    "mgnify_v6_lsu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-lsu/silva-lsu_138.1.tar.gz",
    "mgnify_v6_ssu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-ssu/silva-ssu_138.1.tar.gz",
    "mgnify_v6_its_unite": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/unite/unite_9.0.tar.gz",
    "mgnify_v6_its_itsonedb": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/itsonedb/itsonedb_1.141.tar.gz",
    "mgnify_v6_pr2": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pr2/pr2_5.0.0.tar.gz",
    "test_lsu": "https://zenodo.org/record/8205348/files/test_lsu.tar.gz",
}

DB_names = {
    "mgnify_v5_lsu": "MGnify LSU (v5.0.7) - silva_lsu-20200130",
    "mgnify_v5_ssu": "MGnify SSU (v5.0.7) - silva_ssu-20200130",
    "mgnify_v5_its_unite": "MGnify ITS UNITE (v5.0.7) - UNITE-20200214",
    "mgnify_v5_its_itsonedb": "MGnify ITS ITSonedb (v5.0.7) - ITSoneDB-20200214",
    "mgnify_v6_lsu": "MGnify LSU (v6.0) - silva_lsu-20240702",
    "mgnify_v6_ssu": "MGnify SSU (v6.0) - silva_ssu-20240701",
    "mgnify_v6_its_unite": "MGnify ITS UNITE (v6.0) - UNITE-20240702",
    "mgnify_v6_its_itsonedb": "MGnify ITS ITSonedb (v6.0) - ITSoneDB-20240702",
    "mgnify_v6_pr2": "MGnify PR2 (v6.0) - PR2-20240702",
    "test_lsu": "Trimmed LSU Test DB",
}

def _empty_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)
        return
    for name in os.listdir(path):
        p = os.path.join(path, name)
        if os.path.isdir(p) and not os.path.islink(p):
            shutil.rmtree(p)
        else:
            try:
                os.remove(p)
            except OSError:
                pass

def _copy_all(src_dir, dest_dir):
    for name in os.listdir(src_dir):
        s = os.path.join(src_dir, name)
        d = os.path.join(dest_dir, name)
        if os.path.isdir(s) and not os.path.islink(s):
            if os.path.exists(d):
                shutil.rmtree(d)
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, dest_dir)

def _safe_members(members):
    for m in members:
        mpath = m.name
        norm = mpath.replace("\\", "/")
        if os.path.isabs(mpath) or ".." in norm.split("/"):
            continue
        yield m

def materialize_db_into(output_dir, url):
    # Work entirely under /tmp to avoid job-dir races
    tmp_root = tempfile.mkdtemp(prefix="mapseq_dm_", dir="/tmp")
    temp_extract = tempfile.mkdtemp(prefix="extract_", dir=tmp_root)

    tar_path = wget.download(url, out=tmp_root)
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(temp_extract, members=_safe_members(tar))

    _empty_dir(output_dir)

    entries = os.listdir(temp_extract)
    if len(entries) == 1 and os.path.isdir(os.path.join(temp_extract, entries[0])):
        inner = os.path.join(temp_extract, entries[0])
        vf = os.path.join(inner, "VERSION.txt")
        if os.path.exists(vf):
            try:
                os.remove(vf)
            except OSError:
                pass
        _copy_all(inner, output_dir)
    else:
        vf = os.path.join(temp_extract, "VERSION.txt")
        if os.path.exists(vf):
            try:
                os.remove(vf)
            except OSError:
                pass
        _copy_all(temp_extract, output_dir)

    try:
        shutil.rmtree(tmp_root)
    except Exception:
        pass

def main():
    parser = argparse.ArgumentParser(description="Fetch and register MAPseq DB")
    parser.add_argument("--out", dest="output", required=True, help="Galaxy params/DM JSON path")
    parser.add_argument("--version", dest="version", required=False, help="DB version label (e.g., 6.0)")
    parser.add_argument("--database-type", dest="db_type", required=False, help="DB key (e.g., mgnify_v6_lsu)")
    parser.add_argument("--test", action="store_true", help="Use small test DB")
    args = parser.parse_args()

    # Read Galaxy params to get extra_files_path
    with open(args.output) as fh:
        params = json.load(fh)

    output_dir = params["output_data"][0]["extra_files_path"]
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Resolve db_type/version: prefer CLI, else try params.select_version
    db_type = args.db_type
    version = args.version
    sel = params.get("select_version") or {}
    if not db_type:
        db_type = sel.get("database_type")
    if not version:
        version = sel.get("version")

    if not db_type:
        raise RuntimeError("Missing --database-type and no params['select_version']['database_type'] provided.")
    if db_type not in DB_paths:
        raise KeyError("Unknown database type: {}. Valid keys: {}".format(db_type, ", ".join(sorted(DB_paths.keys()))))

    url = DB_paths["test_lsu"] if args.test else DB_paths[db_type]
    materialize_db_into(output_dir, url)

    # Build DM JSON (value/name pattern you used earlier)
    date_str = datetime.utcnow().strftime("%Y-%m-%d")
    db_value = "{}_from_{}".format(db_type, date_str)
    db_name = DB_names.get(db_type, db_type)
    entry_name = "{} downloaded at {}".format(db_name, date_str)

    data_manager_entry = {
        "data_tables": {
            "mapseq_db": [
                {
                    "value": db_value,
                    "name": entry_name,
                    "path": output_dir,
                    "version": version,
                }
            ]
        }
    }

    with open(args.output, "w") as fh:
        json.dump(data_manager_entry, fh, indent=2, sort_keys=True)

    try:
        count = len(os.listdir(output_dir))
        print("[INFO] Wrote {} items into {}".format(count, output_dir))
    except Exception:
        pass

    print("[SUCCESS] Data Manager JSON written to {}".format(args.output))

if __name__ == "__main__":
    main()
