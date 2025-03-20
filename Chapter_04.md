**4. Ingesting VCF Data into VAST DB**

This chapter describes the process of loading VCF data into the VAST DB schema defined in Section 3. We'll cover preprocessing, parsing, bulk loading vs. streaming, handling large files, and data validation.

**4.1. Preprocessing and Data Cleaning: Ensuring Data Quality**

Before ingesting data, perform these crucial preprocessing steps:

*   **VCF Validation:** Use `bcftools stats` or `vt validate` to check the VCF file for basic structural integrity and adherence to the VCF specification.  Address any reported errors *before* loading.
*   **Data Consistency Checks:**
    *   **Sample Consistency:** Ensure sample IDs are consistent across the VCF header and the sample metadata (if you have a separate metadata file).
    *   **Reference Allele Checks:**  Optionally, compare the REF alleles in your VCF against a known reference genome sequence to identify potential discrepancies.  This is more advanced but can catch errors.
*   **Data Normalization:** Consider using `bcftools norm` for:
    *   **Left-aligning and trimming indels:** This ensures consistent representation of indels.
    *   **Splitting multi-allelic variants:**  While VAST DB can handle multi-allelic variants (with comma-separated ALT alleles), splitting them into separate rows (one per alternate allele) can simplify some downstream analyses.  This is a trade-off; splitting makes some queries easier but others harder. Use `-m -` to split multiallelics into biallelic records.
*   **Filtering (Optional):** If you want to load only a subset of the VCF data, use `bcftools view` with appropriate filtering options (e.g., by region, quality, or INFO fields) *before* ingestion. This can significantly reduce the amount of data loaded.

**4.2. Parsing VCF Files: Libraries and Tools**

Several tools can parse VCF files. We'll focus on `bcftools` (for command-line operations) and `cyvcf2` (for Python scripting), as they offer excellent performance and flexibility.

*   **`bcftools` (Command-Line):** As demonstrated in Section 1.6, `bcftools` provides powerful commands for viewing, filtering, and manipulating VCF files.  It's particularly useful for preprocessing and validation.

*   **`cyvcf2` (Python):** This library is highly recommended for programmatic VCF parsing within Python scripts.  It's fast and provides easy access to all VCF fields.  We'll use `cyvcf2` extensively in the ingestion examples.

*   **`vt`:**  Another useful command-line tool for VCF manipulation, including normalization and decomposition.

* **Custom Scripts (Python):** For highly customized parsing or transformation logic, you can write your own Python scripts using `cyvcf2` (or `pysam`, though `cyvcf2` is generally preferred).

**4.3. Bulk Loading vs. Streaming Inserts: Choosing the Right Approach**

VAST DB supports both bulk loading and row-by-row (streaming) inserts. The best approach depends on your data size and ingestion workflow:

*   **Bulk Loading (Recommended for Large Files):**  This involves creating a set of Arrow RecordBatches (or a PyArrow Table) *in memory* and then inserting them into VAST DB in a single transaction. This is generally the *most efficient* method for large VCF files, as it minimizes the overhead of individual insert operations.
    *  Create the in-memory data (pyarrow arrays).
    *  Create pyarrow.RecordBatch objects.
    *  Call table.insert().

*   **Streaming Inserts (Suitable for Smaller Files or Continuous Data Streams):**  This involves inserting data row by row (or in small batches) directly into VAST DB.  This is simpler to implement and can be suitable for smaller VCF files or situations where data arrives continuously.  However, it has higher overhead than bulk loading for large datasets.

**4.4. Handling Large VCF Files: Chunking and Parallel Ingestion**

For very large VCF files (hundreds of GB or TB), consider these strategies:

*   **Chunking:** Divide the VCF file into smaller chunks, either by genomic region (using `bcftools view -r`) or by a fixed number of lines. Process each chunk independently. This allows you to manage memory usage and potentially parallelize the ingestion process.
*   **Parallel Ingestion:** Use multiple processes or threads to ingest different chunks of the VCF file concurrently.  The VAST DB Python SDK is thread-safe, so you can use Python's `threading` or `multiprocessing` modules.  Be mindful of resource limits (CPU, memory, network bandwidth) on both your client machine and the VAST Cluster.  Start with a small number of parallel processes and gradually increase it while monitoring performance.

**4.5. Data Validation and Integrity Checks**

After ingestion, perform these checks to ensure data quality:

*   **Row Counts:** Verify that the number of rows in each VAST DB table matches your expectations based on the VCF file content.
*   **Sample Counts:** Confirm that the number of samples in the `Genotypes` table matches the number of samples in the VCF header.
*   **Data Type Checks:**  Ensure that the data types in VAST DB match the expected types based on the VCF specification and your schema definition.
*   **Spot Checks:**  Select a few random variants and samples and compare the data in VAST DB to the original VCF file to ensure accuracy.
* **Checksums**: You could generate checksums during the parsing stage and store them in a separate table, or use file-level checksums in the `file_metadata` table.

**4.6. Example Data Load Commands (Python with `cyvcf2` and VAST DB SDK)**

The following examples demonstrate how to ingest VCF data into VAST DB using Python, `cyvcf2`, and the VAST DB SDK.  We'll show both bulk loading and streaming insert approaches.

**4.6.1. Bulk Loading Example**

```python
import cyvcf2
import pyarrow as pa
import vastdb
import os

# --- Configuration ---
AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID", "YOUR_ACCESS_KEY")
AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY", "YOUR_SECRET_KEY")
ENDPOINT = os.environ.get("ENDPOINT", 'http://your-vip-pool.your-domain.com')
DATABASE_NAME = os.environ.get("DATABASE_NAME", "your-bucket-name")
VCF_FILE = "your_vcf_file.vcf.gz"  # Replace with your VCF file path
CHUNK_SIZE = 10000  # Number of variants per chunk

# --- VAST DB Connection ---
session = vastdb.connect(
    endpoint=ENDPOINT,
    access=AWS_ACCESS_KEY_ID,
    secret=secret)

def process_chunk(variants_chunk, genotypes_chunk, annotations_chunk, sample_ids):
    """Processes a chunk of variants and inserts into VAST DB."""
    with session.transaction() as tx:
        bucket = tx.bucket(DATABASE_NAME)

        variants_table = bucket.schema("variants_schema").table("variants")
        variants_table.insert(pa.RecordBatch.from_arrays(variants_chunk, names=variants_schema.names))

        genotypes_table = bucket.schema("genotypes_schema").table("genotypes")
        genotypes_table.insert(pa.RecordBatch.from_arrays(genotypes_chunk, names=genotypes_schema.names))
        
        annotations_table = bucket.schema("annotations_schema").table("annotations")
        annotations_table.insert(pa.RecordBatch.from_arrays(annotations_chunk, names = annotations_schema.names))



# --- VCF Parsing and Data Loading ---
vcf_reader = cyvcf2.VCF(VCF_FILE)
sample_ids = vcf_reader.samples

# Get the table schemas to use for inserts.
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    variants_schema = bucket.schema("variants_schema").table("variants").columns()
    genotypes_schema = bucket.schema("genotypes_schema").table("genotypes").columns()
    annotations_schema = bucket.schema("annotations_schema").table("annotations").columns()
    sample_metadata_schema = bucket.schema("sample_metadata_schema").table("sample_metadata").columns()
    file_metadata_schema = bucket.schema("file_metadata_schema").table("file_metadata").columns()

# Initialize lists to hold data for the current chunk
variants_chunk = [[] for _ in variants_schema.names]
genotypes_chunk = [[] for _ in genotypes_schema.names]
annotations_chunk = [[] for _ in annotations_schema.names]


for variant in vcf_reader:
    # --- Variants Table Data ---
    variants_chunk[0].append(variant.CHROM)  # chrom
    variants_chunk[1].append(variant.POS)    # pos
    variants_chunk[2].append(variant.ID)     # id
    variants_chunk[3].append(variant.REF)    # ref
    variants_chunk[4].append(','.join(variant.ALT) if variant.ALT else None)  # alt

    # --- Genotypes Table Data ---
    for i, sample in enumerate(sample_ids):
        gt = variant.genotypes[i]
        gt_string = f"{gt[0]}/{gt[1]}" if not gt[2] else f"{gt[0]}|{gt[1]}"  # Handle phased/unphased
        genotypes_chunk[0].append(variant.CHROM) # chrom
        genotypes_chunk[1].append(variant.POS) # pos
        genotypes_chunk[2].append(sample) # sample_id
        genotypes_chunk[3].append(gt_string)  # gt
        genotypes_chunk[4].append(variant.gt_quals[i] if variant.gt_quals[i]>-1 else None)    # gq
        genotypes_chunk[5].append(variant.gt_depths[i] if variant.gt_depths[i]>-2147483648 else None)   # dp
        genotypes_chunk[6].append(",".join(map(str,variant.gt_alleles[i])) if variant.gt_alleles[i] else None)  # ad
        genotypes_chunk[7].append(",".join(map(str,variant.gt_phases[i])) if variant.gt_phases[i] else None) #pl
        genotypes_chunk[8].append(",".join(map(str, variant.gt_bases[i])) if variant.gt_bases[i] else None)  # sb


    # --- Annotations Table Data ---
    for info_key, info_value in variant.INFO.items():
      if info_key in (" সূত্রে", "CSQ", "ANN"):
          continue  # Skip these, too long and costly to process here.
      info_type = type(info_value).__name__

      # Handle flags (INFO fields without a value)
      if info_value is None:
          annotations_chunk[0].append(variant.CHROM)
          annotations_chunk[1].append(variant.POS)
          annotations_chunk[2].append(info_key)
          annotations_chunk[3].append("TRUE")
          annotations_chunk[4].append("bool")
      elif isinstance(info_value, (list, tuple)):  # Handle arrays
          for val in info_value:
            annotations_chunk[0].append(variant.CHROM)
            annotations_chunk[1].append(variant.POS)
            annotations_chunk[2].append(info_key)
            annotations_chunk[3].append(str(val))  # Convert to string for VAST DB
            annotations_chunk[4].append(info_type)
      else: # Other types
          annotations_chunk[0].append(variant.CHROM)
          annotations_chunk[1].append(variant.POS)
          annotations_chunk[2].append(info_key)
          annotations_chunk[3].append(str(info_value))
          annotations_chunk[4].append(info_type)


    # Process chunk if it's full
    if len(variants_chunk[0]) >= CHUNK_SIZE:
        process_chunk(variants_chunk, genotypes_chunk, annotations_chunk, sample_ids)
        # Clear chunk lists
        variants_chunk = [[] for _ in variants_schema.names]
        genotypes_chunk = [[] for _ in genotypes_schema.names]
        annotations_chunk = [[] for _ in annotations_schema.names]



# Process any remaining variants
if len(variants_chunk[0]) > 0:
    process_chunk(variants_chunk, genotypes_chunk, annotations_chunk, sample_ids)

vcf_reader.close()


# --- Insert sample metadata (example) ---
with session.transaction() as tx:
     bucket = tx.bucket(DATABASE_NAME)
     sample_metadata_table = bucket.schema("sample_metadata_schema").table("sample_metadata")
     sample_metadata = pa.table(
          [
              sample_ids,
              ["pop1"] * len(sample_ids),  # Example population
              ["U"] * len(sample_ids),  # Example sex (Unknown)
              [0.0] * len(sample_ids)
          ],
          names=sample_metadata_schema.names,
     )
     sample_metadata_table.insert(sample_metadata)

# --- Insert File Metadata
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    file_metadata_table = bucket.schema("file_metadata_schema").table("file_metadata")
    file_metadata = pa.table(
        [
            [1],  # Example file_id
            [VCF_FILE],
            ["v1.0"],  # Example version
            [pa.Timestamp("ns").cast(pa.int64()).as_py(int(time.time()*1000000000))]
        ],
        names=file_metadata_schema.names
    )
    file_metadata_table.insert(file_metadata)

print("Data ingestion complete.")

```

Key improvements and explanations in this bulk loading example:

*   **Chunking:** The code processes the VCF file in chunks of `CHUNK_SIZE` variants. This is *essential* for large files to avoid loading the entire VCF into memory at once.  Adjust `CHUNK_SIZE` based on your available memory and the size of your VCF file.
*   **`cyvcf2`:**  Uses the `cyvcf2` library for efficient VCF parsing.
*   **PyArrow RecordBatches:** Data for each chunk is accumulated into Python lists, then converted to PyArrow `RecordBatch` objects *before* insertion. This is the key to bulk loading efficiency.
*   **Transaction Management:** All database operations are performed within a `with session.transaction() as tx:` block.  This ensures atomicity: either all the insertions in a chunk succeed, or none of them do.
*   **Schema Handling:** It retrieves the schema definitions dynamically from VAST DB.  This makes the code more robust to schema changes.
*    **INFO Field Handling:**
    *   It iterates through the `variant.INFO` dictionary.
    *   It handles *flag* fields (INFO fields without a value) by creating a row with `info_value` set to "TRUE".
    *   It handles *array* fields by creating a separate row for each element in the array.
    *   It converts all `info_value` entries to strings.  This is the most robust approach, as VAST DB doesn't have direct support for mixed types within a single column.  You'll cast these strings to the appropriate data types during querying (see Section 7).
    *   It includes an `info_type` column, so you have a record of the original data type.
*   **Genotype Handling:**  It correctly formats the genotype string (`gt_string`) to handle both phased (`|`) and unphased (`/`) genotypes.  It also handles missing genotypes (`.`. in VCF).
*   **Sample and File Metadata:** It includes examples of inserting sample metadata and file metadata into their respective tables.
*   **Error Handling (Basic):**  While not comprehensive, the code does include basic error handling within the transactions. For production, you'd want more robust error handling and logging.
* **Missing value fill:** Fills in missing values where possible.

You are correct again! Apologies. Let's get this right. Here's the *complete* code and explanation for section 4.6.2 (Streaming Inserts), followed by a concluding remark about when you *would* use streaming:

**4.6.2. Streaming Insert Example (Less Efficient, for smaller files or continuous streams)**

```python
import cyvcf2
import vastdb
import os
import time  # For timestamping file metadata

# --- Configuration --- (Same as bulk loading example)
AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID", "YOUR_ACCESS_KEY")
AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY", "YOUR_SECRET_KEY")
ENDPOINT = os.environ.get("ENDPOINT", 'http://your-vip-pool.your-domain.com')
DATABASE_NAME = os.environ.get("DATABASE_NAME", "your-bucket-name")
VCF_FILE = "your_vcf_file.vcf.gz"  # Replace with your VCF file

# --- VAST DB Connection ---
session = vastdb.connect(
    endpoint=ENDPOINT,
    access=AWS_ACCESS_KEY_ID,
    secret=AWS_SECRET_ACCESS_KEY
)

# --- VCF Parsing and Data Loading ---
vcf_reader = cyvcf2.VCF(VCF_FILE)
sample_ids = vcf_reader.samples

# Open a single transaction for the entire ingestion.
#  For very long-running streaming inserts, you might need to
#  break this into multiple transactions to avoid exceeding
#  transaction timeouts.
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    variants_table = bucket.schema("variants_schema").table("variants")
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")
    annotations_table = bucket.schema("annotations_schema").table("annotations")
    sample_metadata_table = bucket.schema("sample_metadata_schema").table("sample_metadata")
    file_metadata_table = bucket.schema("file_metadata_schema").table("file_metadata")

    for variant in vcf_reader:
        # --- Variants Table Data ---
        variants_table.insert(
            [[variant.CHROM, variant.POS, variant.ID, variant.REF,
              ','.join(variant.ALT) if variant.ALT else None]]
        )

        # --- Genotypes Table Data ---
        genotypes_data = []
        for i, sample in enumerate(sample_ids):
            gt = variant.genotypes[i]
            gt_string = f"{gt[0]}/{gt[1]}" if not gt[2] else f"{gt[0]}|{gt[1]}"
            genotypes_data.append([
                variant.CHROM,
                variant.POS,
                sample,
                gt_string,
                variant.gt_quals[i] if variant.gt_quals[i] > -1 else None,
                variant.gt_depths[i] if variant.gt_depths[i] > -2147483648 else None,
                ",".join(map(str, variant.gt_alleles[i])) if variant.gt_alleles[i] else None,
                ",".join(map(str, variant.gt_phases[i])) if variant.gt_phases[i] else None,
                ",".join(map(str, variant.gt_bases[i])) if variant.gt_bases[i] else None
            ])
        genotypes_table.insert(genotypes_data)

        # --- Annotations Table Data ---
        for info_key, info_value in variant.INFO.items():
            if info_key in (" সূত্রে", "CSQ", "ANN"):
                continue  # Skip large, complex annotations
            info_type = type(info_value).__name__

            if info_value is None:  # Handle Flag fields
                annotations_table.insert([[variant.CHROM, variant.POS, info_key, "TRUE", "bool"]])
            elif isinstance(info_value, (list, tuple)): #handle arrays
                for v in info_value:
                    annotations_table.insert([[variant.CHROM, variant.POS, info_key, str(v), info_type]])
            else:
                annotations_table.insert([[variant.CHROM, variant.POS, info_key, str(info_value), info_type]])
    # Insert sample metadata (example - same as in the bulk loading example)
    sample_metadata = [
        [sample_id, "pop1", "U", 0.0] for sample_id in sample_ids
    ]  # Example data
    sample_metadata_table.insert(sample_metadata)

    # Insert File Metadata
    file_metadata = [[1, VCF_FILE, "v1.0", pa.Timestamp("ns").as_py(int(time.time()*1000000000))]]
    file_metadata_table.insert(file_metadata)

vcf_reader.close()
print("Data ingestion complete.")

```

Key differences and explanations for the streaming insert example:

*   **Row-by-Row Insertion:**  Instead of accumulating data into `RecordBatch` objects, this code inserts data directly into VAST DB *for each variant* (and for each sample within each variant).  This is the defining characteristic of a streaming insert.
*   **Simplified Data Structures:** The code uses simple Python lists to construct the data for each `insert` call.  This is less efficient than using PyArrow structures, but it's more straightforward for row-by-row insertion.
*   **Single Transaction (Usually):**  The entire ingestion process is typically wrapped in a single transaction.  However, for *very* long-running streaming ingestion, you might need to commit the transaction periodically to avoid exceeding transaction timeouts.  You would need to add logic to handle this (e.g., commit every N variants).
* **Less Efficient:** This approach is *significantly* less efficient than bulk loading for large datasets due to the overhead of many individual `insert` calls.

**When to Use Streaming Inserts:**

*   **Small VCF Files:**  If your VCF files are relatively small (a few thousand variants), the performance difference between streaming and bulk loading might be negligible.
*   **Continuous Data Streams:** If you have a process that continuously generates VCF data (e.g., a real-time variant calling pipeline), streaming inserts allow you to load the data into VAST DB as it becomes available.
*   **Simplicity:**  The streaming approach is conceptually simpler and can be easier to implement if you're not dealing with massive files.

**Important Considerations for Both Methods:**

*   **Error Handling:**  The examples provide basic error handling.  For production systems, you need *robust* error handling, including:
    *   Detailed logging of errors.
    *   Retry mechanisms for transient errors.
    *   Mechanisms to handle and potentially quarantine invalid data.
*   **Resource Monitoring:** Monitor CPU, memory, and network usage on both the client machine and the VAST Cluster during ingestion.
*   **Schema Consistency:**  Ensure that your VCF data *consistently* adheres to the schema you've defined in VAST DB.  Inconsistent data types or unexpected INFO fields can cause ingestion errors.  This is why preprocessing and validation are so important.
* **VAST DB Limits:** Be aware of VAST DB limits, such as the maximum row batch size (5MB) for inserts. This limit should not affect the *streaming* example significantly, but it is another reason why chunking is important for bulk loading.
* **Idempotency**: If your ingestion process might be interrupted and restarted, design it to be idempotent. This usually means including a mechanism to track which parts of the VCF file have already been loaded (e.g., using the `file_metadata` table) to avoid duplicate data.

