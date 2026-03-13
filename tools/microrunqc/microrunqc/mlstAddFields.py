#!/usr/bin/env 

import sys
import csv

def find_index(headers, term):
    try:
        return headers.index(term)
    except ValueError:
        return -1

def main(mlst_file, db_path=None):
    with open(mlst_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        mlstout = next(reader)

    schema = mlstout[1]
    mlstST = mlstout[2]

    # Return the output without appending if schema equals "-"
    if schema == "-":
        print("\t".join(mlstout))
        return 
        
    if db_path is None:
        # If no database path is provided, find it using an external command
        # This requires the 'mlst' command to be installed and available in the path
        import subprocess
        mlstdesc = subprocess.check_output(['mlst', '-h']).decode()
        db_pubmlst = [line for line in mlstdesc.split('\n') if 'db/pubmlst' in line]
        if db_pubmlst:
            mlstloc = db_pubmlst[0].split("'")[1].replace("bin/..", "")
        else:
            raise Exception("Could not find MLST database location.")
    else:
        mlstloc = db_path

    mlst_file_path = f"{mlstloc}/{schema}/{schema}.txt"

    schema_dict = {}
    with open(mlst_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)

        clonal = find_index(headers, 'clonal_complex')
        cc = find_index(headers, 'CC')
        lineage = find_index(headers, 'Lineage')
        species = find_index(headers, 'species')

        for line in reader:
            desc = []
            if clonal > -1 and line[clonal]:
                desc.append(f"clonal_complex={line[clonal]}")
            if cc > -1 and line[cc]:
                desc.append(f"CC={line[cc]}")
            if lineage > -1 and line[lineage]:
                desc.append(f"Lineage={line[lineage]}")
            if species > -1 and line[species]:
                desc.append(f"species={line[species]}")
            schema_dict[line[0]] = ','.join(desc)

    output = mlstout[:3]
    if mlstST in schema_dict:
        output.append(schema_dict[mlstST])
    else:
        output.append("-")
    output.extend(mlstout[3:])

    print("\t".join(output))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python mlstAddFields.py <mlst_file> [db_path]")
        sys.exit(1)

    mlst_file = sys.argv[1]
    db_path = sys.argv[2] if len(sys.argv) > 2 else None

    main(mlst_file, db_path)

