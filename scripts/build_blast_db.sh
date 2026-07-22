#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ZIP_DIR="${ROOT_DIR}/FLU_DB"
ZIP="${ZIP_DIR}/fluA_reference.fasta.zip"
DB_DIR="${ZIP_DIR}/fluA_db"
FASTA="${DB_DIR}/fluA_reference.fasta"
DB_PREFIX="${DB_DIR}/fluA_db"
TMP_DIR="${DB_DIR}/_extract_tmp"

GITHUB_REPO="${GITHUB_REPO:-ZooPhy/apgap-influenza-pipeline}"
RELEASE_TAG="${RELEASE_TAG:-v0.1.0}"
ASSET_NAME="${ASSET_NAME:-fluA_reference.fasta.zip}"
ZIP_URL="${ZIP_URL:-https://github.com/${GITHUB_REPO}/releases/download/${RELEASE_TAG}/${ASSET_NAME}}"
BLAST_IMAGE="${BLAST_IMAGE:-ncbi/blast}"

mkdir -p "${ZIP_DIR}" "${DB_DIR}"

if [[ ! -f "${ZIP}" ]]; then
    echo "Compressed reference not found locally: ${ZIP}"
    echo "Attempting to download release asset: ${ASSET_NAME}"
    echo "Repository: ${GITHUB_REPO}"
    echo "Release tag: ${RELEASE_TAG}"

    download_ok=false

    # Prefer GitHub CLI when available. This supports private repositories
    # when the user has authenticated with: gh auth login
    if command -v gh >/dev/null 2>&1; then
        echo "Trying GitHub CLI download..."
        if gh release download "${RELEASE_TAG}" \
            --repo "${GITHUB_REPO}" \
            --pattern "${ASSET_NAME}" \
            --dir "${ZIP_DIR}" \
            --clobber; then
            if [[ -s "${ZIP}" ]]; then
                download_ok=true
            else
                echo "GitHub CLI completed, but the expected asset was not found:"
                echo "  ${ZIP}"
            fi
        else
            rm -f "${ZIP}"
            echo "GitHub CLI download failed; trying direct URL."
        fi
    fi

    # Direct curl works for public repositories and public release assets.
    if [[ "${download_ok}" != true ]]; then
        echo "Downloading from: ${ZIP_URL}"
        if curl -fL --retry 3 --retry-delay 2 -o "${ZIP}" "${ZIP_URL}"; then
            download_ok=true
        else
            rm -f "${ZIP}"
        fi
    fi

    if [[ "${download_ok}" != true ]]; then
        cat >&2 <<ERROR
ERROR: Unable to download ${ASSET_NAME}.

The release URL is valid only when all of the following are true:
  1. The repository ${GITHUB_REPO} is publicly accessible, or GitHub CLI is authenticated.
  2. A release with tag ${RELEASE_TAG} exists.
  3. That release contains an asset named exactly ${ASSET_NAME}.

Fix one of the following:
  - Create/publish the ${RELEASE_TAG} release and attach ${ASSET_NAME}.
  - For a private repository, run 'gh auth login' and rerun this script.
  - Place the archive manually at:
      ${ZIP}
  - Override the source URL:
      ZIP_URL='https://example.org/${ASSET_NAME}' ./scripts/build_blast_db.sh
  - Override the release tag:
      RELEASE_TAG='vX.Y.Z' ./scripts/build_blast_db.sh
ERROR
        exit 1
    fi
fi

if [[ ! -s "${ZIP}" ]]; then
    echo "ERROR: Downloaded archive is missing or empty: ${ZIP}" >&2
    exit 1
fi

echo "Extracting influenza reference..."
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"
unzip -oq "${ZIP}" -d "${TMP_DIR}"

FASTA_SRC="$(find "${TMP_DIR}" -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \) | head -n 1 || true)"

if [[ -z "${FASTA_SRC}" ]]; then
    echo "ERROR: No FASTA file found inside ${ZIP}" >&2
    rm -rf "${TMP_DIR}"
    exit 1
fi

cp -f "${FASTA_SRC}" "${FASTA}"

# Remove old database files before rebuilding so version 4 and version 5
# database components cannot be mixed.
find "${DB_DIR}" -maxdepth 1 -type f \
    \( -name 'fluA_db.n*' -o -name 'fluA_db.*db' \) \
    -delete

if command -v docker >/dev/null 2>&1; then
    echo "Building BLAST database version 4 using Docker..."

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
    echo "Building BLAST database version 4 using Apptainer..."

    apptainer exec \
        --bind "${ROOT_DIR}:/data" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/data/FLU_DB/fluA_db/fluA_reference.fasta" \
            -dbtype nucl \
            -blastdb_version 4 \
            -out "/data/FLU_DB/fluA_db/fluA_db"

elif command -v singularity >/dev/null 2>&1; then
    echo "Building BLAST database version 4 using Singularity..."

    singularity exec \
        --bind "${ROOT_DIR}:/data" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/data/FLU_DB/fluA_db/fluA_reference.fasta" \
            -dbtype nucl \
            -blastdb_version 4 \
            -out "/data/FLU_DB/fluA_db/fluA_db"

else
    echo "ERROR: Neither Docker, Apptainer, nor Singularity was found." >&2
    echo "Install one of these container runtimes before running this script." >&2
    rm -rf "${TMP_DIR}"
    exit 1
fi

rm -rf "${TMP_DIR}"

echo
echo "BLAST database successfully created:"
echo "  ${DB_PREFIX}"
