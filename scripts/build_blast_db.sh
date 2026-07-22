#!/usr/bin/env bash
set -euo pipefail

# Directory containing this script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Repository root
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ZIP_DIR="${ROOT_DIR}/FLU_DB"
ZIP="${ZIP_DIR}/fluA_reference.fasta.zip"
DB_DIR="${ZIP_DIR}/fluA_db"
FASTA="${DB_DIR}/fluA_reference.fasta"
DB_PREFIX="${DB_DIR}/fluA_db"
TMP_DIR="${DB_DIR}/_extract_tmp"

# Update this when you publish the release asset
RELEASE_TAG="${RELEASE_TAG:-v0.1.0}"
ZIP_URL="${ZIP_URL:-https://github.com/ZooPhy/apgap-influenza-pipeline/releases/download/${RELEASE_TAG}/fluA_reference.fasta.zip}"

# BLAST container image can be overridden if needed
BLAST_IMAGE="${BLAST_IMAGE:-ncbi/blast}"

mkdir -p "${DB_DIR}"

if [[ ! -f "${ZIP}" ]]; then
    echo "Compressed reference not found locally."
    echo "Downloading from: ${ZIP_URL}"
    curl -fL --retry 3 --retry-delay 2 -o "${ZIP}" "${ZIP_URL}"
fi

if [[ ! -f "${ZIP}" ]]; then
    echo "ERROR: ${ZIP} not found."
    exit 1
fi

echo "Extracting influenza reference..."
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

unzip -oq "${ZIP}" -d "${TMP_DIR}"

FASTA_SRC="$(find "${TMP_DIR}" -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \) | head -n 1 || true)"

if [[ -z "${FASTA_SRC}" ]]; then
    echo "ERROR: No FASTA file found inside ${ZIP}"
    rm -rf "${TMP_DIR}"
    exit 1
fi

cp -f "${FASTA_SRC}" "${FASTA}"

if command -v docker >/dev/null 2>&1; then
    echo "Building BLAST database using Docker..."

    docker run --rm --platform linux/amd64 \
        -v "${ROOT_DIR}:/data" \
        -w /data \
        "${BLAST_IMAGE}" \
        makeblastdb \
            -in "/data/FLU_DB/fluA_db/fluA_reference.fasta" \
            -dbtype nucl \
            -blastdb_version 4 \
            -out "/data/FLU_DB/fluA_db/fluA_db"

elif command -v apptainer >/dev/null 2>&1; then
    echo "Building BLAST database using Apptainer..."

    apptainer exec \
        --bind "${ROOT_DIR}:/data" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/data/FLU_DB/fluA_db/fluA_reference.fasta" \
            -dbtype nucl \
            -out "/data/FLU_DB/fluA_db/fluA_db"

elif command -v singularity >/dev/null 2>&1; then
    echo "Building BLAST database using Singularity..."

    singularity exec \
        --bind "${ROOT_DIR}:/data" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/data/FLU_DB/fluA_db/fluA_reference.fasta" \
            -dbtype nucl \
            -out "/data/FLU_DB/fluA_db/fluA_db"

else
    echo "ERROR: Neither Docker, Apptainer, nor Singularity was found."
    echo "Install one of these container runtimes before running this script."
    rm -rf "${TMP_DIR}"
    exit 1
fi

rm -rf "${TMP_DIR}"

echo
echo "BLAST database successfully created:"
echo "  ${DB_PREFIX}"