#!/bin/bash

# Downloads necessary data to ingest into the genewalkR DB

if [ "$#" -ne 2 ] || [ "$1" != "-t" ]; then
    echo "Usage: $0 -t <target_folder>"
    exit 1
fi

TARGET_FOLDER="$2"
OT_URL="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"

mkdir -p "$TARGET_FOLDER/targets"
mkdir -p "$TARGET_FOLDER/reactome"
mkdir -p "$TARGET_FOLDER/interactions"
mkdir -p "$TARGET_FOLDER/interaction_pairs"
mkdir -p "$TARGET_FOLDER/obo_files"

download_parquet_files() {
    local url="$1"
    local dest="$2"
    
    echo "Downloading from $url"
    local files=$(curl -s "$url" | grep -oE 'href="[^"]+.parquet"' | cut -d'"' -f2)
    
    for file in $files; do
        echo "Downloading $file..."
        wget "${url}${file}" -P "$dest" --no-verbose
    done
}

download_parquet_files "$OT_URL/target/" "$TARGET_FOLDER/targets"

download_parquet_files "$OT_URL/reactome/" "$TARGET_FOLDER/reactome"

download_parquet_files "$OT_URL/interaction_evidence/" "$TARGET_FOLDER/interactions"

# since 26.09 the participants live here, joined on interactionId
download_parquet_files "$OT_URL/interaction/" "$TARGET_FOLDER/interaction_pairs"

echo "Downloading GO obo file..."
wget "https://purl.obolibrary.org/obo/go.obo" -P "$TARGET_FOLDER/obo_files" --no-verbose

echo "Done."