"""Fetch NCBI SRA runinfo rows with E-utilities efetch."""
import csv
import http.client
import time
import urllib.error
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REQUESTS_PER_SECOND = {False: 3.0, True: 10.0}
# efetch sometimes answers 400 to a valid request. The same request then succeeds.
RETRY_STATUS = {400, 429, 500, 502, 503, 504}
MAX_BACKOFF_SECONDS = 60

# efetch sends no header row, so the column names are fixed here and parse() checks the width of each row.
COLUMNS = [
    "Run", "ReleaseDate", "LoadDate", "spots", "bases", "spots_with_mates",
    "avgLength", "size_MB", "AssemblyName", "download_path", "Experiment",
    "LibraryName", "LibraryStrategy", "LibrarySelection", "LibrarySource",
    "LibraryLayout", "InsertSize", "InsertDev", "Platform", "Model", "SRAStudy",
    "BioProject", "Study_Pubmed_id", "ProjectID", "Sample", "BioSample",
    "SampleType", "TaxID", "ScientificName", "SampleName", "g1k_pop_code",
    "source", "g1k_analysis_group", "Subject_ID", "Sex", "Disease", "Tumor",
    "Affection_Status", "Analyte_Type", "Histological_Type", "Body_Site",
    "CenterName", "Submission", "dbgap_study_accession", "Consent", "RunHash",
    "ReadHash",
]


class Pacer:
    """Wait between requests so that they stay below the NCBI limit."""

    def __init__(self, with_api_key):
        self.gap = 1.0 / REQUESTS_PER_SECOND[with_api_key]
        self.last = 0.0

    def wait(self):
        delay = self.gap - (time.monotonic() - self.last)
        if delay > 0:
            time.sleep(delay)
        self.last = time.monotonic()


def fetch(accessions, credentials, pacer, retries=8, timeout=120):
    """Return the runinfo CSV text for a list of run accessions. Raise RuntimeError after the last retry."""
    params = {"db": "sra", "rettype": "runinfo", "retmode": "csv", "id": ",".join(accessions)}
    params.update({key: value for key, value in credentials.items() if value})
    body = urllib.parse.urlencode(params).encode()
    for attempt in range(retries + 1):
        pacer.wait()
        try:
            with urllib.request.urlopen(urllib.request.Request(EFETCH, data=body), timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as error:
            failure = error
            if error.code not in RETRY_STATUS:
                break
        except (urllib.error.URLError, http.client.HTTPException, OSError) as error:
            failure = error
        if attempt < retries:
            time.sleep(min(2 ** attempt, MAX_BACKOFF_SECONDS))
    raise RuntimeError(f"efetch failed: {failure!r}")


def parse(text):
    """Return the rows of a runinfo response. Raise ValueError if a row does not have len(COLUMNS) fields."""
    rows = [row for row in csv.reader(text.splitlines()) if any(field.strip() for field in row)]
    if rows and rows[0][0] == "Run":
        rows = rows[1:]
    for row in rows:
        if len(row) != len(COLUMNS):
            raise ValueError(f"runinfo row with {len(row)} fields, expected {len(COLUMNS)}: {row[:3]}")
    return rows
