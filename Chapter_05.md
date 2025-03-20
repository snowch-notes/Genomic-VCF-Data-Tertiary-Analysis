Okay, here's the content for Chapter 5, "Performing Tertiary Analysis with VAST DB," covering basic and advanced queries, statistical analysis, and snapshot usage, all within the context of the VCF data schema and using the VAST DB Python SDK:

**5. Performing Tertiary Analysis with VAST DB**

This chapter demonstrates how to query and analyze VCF data stored in VAST DB, leveraging its columnar storage and MPP architecture for efficient tertiary analysis. We'll cover basic data retrieval, advanced filtering and aggregation, statistical analysis techniques, and the use of snapshots.

**5.1. Basic Queries: Retrieving Variants and Genotypes**

These queries demonstrate fundamental data access patterns. We'll use the VAST DB Python SDK and, for clarity, show the equivalent SQL-like concepts where applicable.

**5.1.1. Filtering by Chromosome, Position, and ID**

This is a fundamental query for retrieving variants within a specific genomic region.

```python
import pyarrow as pa
import vastdb
import os
from ibis import _

# --- Connection (same as before) ---
AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID", "YOUR_ACCESS_KEY")
AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY", "YOUR_SECRET_KEY")
ENDPOINT = os.environ.get("ENDPOINT", 'http://your-vip-pool.your-domain.com')
DATABASE_NAME = os.environ.get("DATABASE_NAME", "your-bucket-name")

session = vastdb.connect(
    endpoint=ENDPOINT,
    access=AWS_ACCESS_KEY_ID,
    secret=AWS_SECRET_ACCESS_KEY
)
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    variants_table = bucket.schema("variants_schema").table("variants")

    # Retrieve variants on chromosome 22 between positions 17000000 and 17001000
    reader = variants_table.select(
        predicate=(
            (_.chrom == "chr22") & (_.pos >= 17000000) & (_.pos <= 17001000)
        )
    )
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())

    # Retrieve variant by ID (assuming 'rs123' exists)
    reader = variants_table.select(predicate=_.id == "rs123")
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())

```

*   **Explanation:**
    *   We use `variants_table.select()` to initiate the query.
    *   The `predicate` parameter specifies the filtering conditions using Ibis expressions.  `_.chrom`, `_.pos`, and `_.id` refer to the respective columns in the `variants` table.
    *   `&` represents the logical AND operator.
    *  The result is a `pyarrow.RecordBatchReader`, and we use list comprehension with `pa.Table.from_batches()` to read all resulting batches into a PyArrow Table, which is then converted to a Pandas DataFrame for easy display.  For very large result sets, you would typically process the `RecordBatchReader` in chunks to avoid loading everything into memory at once.

**5.1.2. Selecting Specific Samples and INFO Fields**

This demonstrates how to retrieve genotype data for specific samples and extract specific INFO field values.

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")
    annotations_table = bucket.schema("annotations_schema").table("annotations")

    # Get genotypes for samples 'sample1' and 'sample2' on chromosome 22, position 17000050
    reader = genotypes_table.select(
        columns=['sample_id', 'gt', 'gq', 'dp'],  # Select specific columns
        predicate=(
            (_.chrom == "chr22") & (_.pos == 17000050) &
            _.sample_id.isin(["sample1", "sample2"])
        )
    )
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())

    # Get the 'AF' (allele frequency) INFO field for variants on chromosome 22, position 17000050
    reader = annotations_table.select(
        columns = ['info_key','info_value'],
        predicate=(
            (_.chrom == "chr22") & (_.pos == 17000050) & (_.info_key == "AF")
        )
    )
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())

```

*   **Explanation:**
    *   `columns=['sample_id', 'gt', 'gq']`:  The `columns` parameter specifies which columns to retrieve (projection pushdown).  This is *crucial* for performance, as it avoids transferring unnecessary data.
    *   `_.sample_id.isin(["sample1", "sample2"])`:  The `isin()` operator efficiently filters for rows where `sample_id` is one of the specified values.
    *   The second query demonstrates retrieving a specific INFO field (`AF`) by filtering on `info_key`.

**5.2. Advanced Queries**

These examples demonstrate more complex filtering, aggregation, and joining operations.

**5.2.1. Calculating Allele Frequencies and Counts**

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")

    # Calculate allele frequencies for variants on chromosome 22
    #  (This example shows a simplified calculation; for real-world use,
    #   you'd need to handle missing genotypes and multi-allelic variants
    #   more carefully).
    reader = genotypes_table.select(
        columns=['chrom', 'pos', 'gt'],
        predicate=_.chrom == "chr22"
    )
    genotype_data = pa.Table.from_batches([b for b in reader]).to_pandas()

    # Basic allele count calculation (assuming biallelic and no missing data)
    def calculate_allele_counts(df):
        counts = {'0': 0, '1': 0}
        for gt in df['gt']:
            if gt is not None:  # Handle missing genotypes
                alleles = gt.replace('|', '/').split('/')  # Handle phased/unphased
                for allele in alleles:
                    if allele in counts:
                        counts[allele] += 1
        return counts

    allele_counts = genotype_data.groupby(['chrom', 'pos']).apply(calculate_allele_counts).reset_index(name='counts')
    print(allele_counts)
    # To get the frequencies, divide by the total number of called alleles.


```

*   **Explanation:**
    *   This example demonstrates a *simplified* allele frequency calculation.  It retrieves all genotypes for chromosome 22, then uses Pandas' `groupby()` and a custom function (`calculate_allele_counts`) to count the occurrences of each allele (0 and 1) at each position.
    *   **Important:** This simplified example doesn't handle missing genotypes (`.`. in VCF) or multi-allelic variants correctly. A production-ready implementation would need to account for these.
    *   This approach loads all relevant genotypes into memory. For very large datasets, consider using VAST DB's aggregation capabilities directly (if supported in future versions) or breaking the calculation into smaller chunks.

**5.2.2. Identifying Rare Variants**

```python
# (Requires the allele frequency calculation from 5.2.1)

# Assuming 'allele_counts' DataFrame is available from the previous example

# Filter for variants with a minor allele frequency (MAF) less than 0.01
# (This assumes 'allele_counts' contains '0' and '1' counts and is biallelic)
def calculate_maf(counts):
  total_alleles = 0
  if '0' in counts:
      total_alleles += counts['0']
  if '1' in counts:
    total_alleles += counts['1']
  if total_alleles == 0:
    return 0

  if '0' not in counts:
    return 1.0
  if '1' not in counts:
    return 0.0
  minor_allele_count = min(counts['0'], counts['1'])
  return float(minor_allele_count) / total_alleles


rare_variants = allele_counts[allele_counts['counts'].apply(calculate_maf) < 0.01]
print(rare_variants)

```

*   **Explanation:** This builds on the previous example. After calculating allele counts, it filters the results to identify variants where the minor allele frequency (MAF) is below a threshold (0.01 in this case). The `calculate_maf` function computes the minor allele frequency, handling cases with one or two observed alleles.

**5.2.3. Filtering by Genotype Quality and Depth**

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")

    # Retrieve genotypes with GQ >= 30 and DP >= 10 on chromosome 22
    reader = genotypes_table.select(
        predicate=(
            (_.chrom == "chr22") & (_.gq >= 30) & (_.dp >= 10)
        )
    )
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())
```

*   **Explanation:** This demonstrates filtering based on genotype quality (`gq`) and read depth (`dp`), common quality control steps in genomic analysis.

**5.2.4. Joining with Annotation Data**

This shows how to combine information from the `Genotypes` and `Annotations` tables.

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")
    annotations_table = bucket.schema("annotations_schema").table("annotations")
    
    # Retrieve genotypes and AF (allele frequency) for a specific sample on chr22

    # First, fetch sample genotypes.  This is much faster than joining.
    genotypes_reader = genotypes_table.select(
        columns=['chrom', 'pos', 'gt'],
        predicate=(
            (_.chrom == "chr22") & (_.sample_id == "sample1")
        )
    )
    genotypes = pa.Table.from_batches([b for b in genotypes_reader]).to_pandas()

    # Then, fetch the annotations.
    annotations_reader = annotations_table.select(
        columns=['chrom', 'pos', 'info_key', 'info_value'],
        predicate = (_.chrom == 'chr22') & (_.info_key == "AF")
    )
    annotations = pa.Table.from_batches([b for b in annotations_reader]).to_pandas()

    # Now, merge them in pandas.
    merged_data = genotypes.merge(annotations, on=['chrom', 'pos'], how='left')
    print(merged_data)
```

*   **Explanation:** This example shows a common pattern: retrieving data from two tables and combining them.  Because VAST DB's `select` method returns a `RecordBatchReader`, it's very efficient to combine result sets using Pandas' `merge` function (or similar operations in other data analysis libraries).
* **Important:**  While SQL-like `JOIN` operations aren't directly supported in the current VAST DB Python SDK *within* the database query itself, you can efficiently perform joins *after* retrieving the data into Pandas DataFrames (or other data structures). This leverages VAST DB's columnar storage and predicate pushdown for efficient data retrieval and then performs the join in memory, which is often quite fast for the result sets typical of tertiary analysis.

**5.2.5. Performing Region-Based Queries (using genomic intervals)**

This expands on 5.1.1 to show how to query for multiple regions efficiently.

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    variants_table = bucket.schema("variants_schema").table("variants")
    # Retrieve variants within multiple regions on chromosome 22
    regions = [
        (17000000, 17001000),
        (18000000, 18002000),
        (19000000, 19001500),
    ]

    # Build the predicate using Ibis's OR operator
    predicate = None
    for start, end in regions:
        region_predicate = (_.chrom == "chr22") & (_.pos >= start) & (_.pos <= end)
        predicate = region_predicate if predicate is None else predicate | region_predicate
    
    reader = variants_table.select(predicate=predicate)
    result_table = pa.Table.from_batches([b for b in reader])
    print(result_table.to_pandas())
```

*   **Explanation:**  This example shows how to query for variants within *multiple* genomic regions using a single query. The `predicate` is constructed by combining individual region predicates using the Ibis `|` (OR) operator. This is much more efficient than issuing separate queries for each region.

**5.3. Statistical Analysis**

VAST DB, in its current form, does not have built-in functions for complex statistical calculations like association tests or LD calculation.  However, its fast data retrieval makes it an excellent *data source* for statistical analysis performed in external tools like R or Python (with libraries like `scipy.stats`, `statsmodels`, or specialized genomic analysis packages).

**5.3.1. Implementing Association Tests (e.g., Chi-squared, Fisher's Exact)**

```python
import scipy.stats as stats
import pandas as pd

with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")
    # Example: Chi-squared test for a single variant on chromosome 22, position 17000050
    #  (Assuming you have case/control status in your sample metadata)

    # 1. Retrieve genotype data for the variant
    reader = genotypes_table.select(
       columns=['sample_id', 'gt'],
       predicate=(
           (_.chrom == "chr22") & (_.pos == 17000050)
       )
    )
    genotype_data = pa.Table.from_batches([b for b in reader]).to_pandas()


    # 2. Merge with sample metadata (assuming you have a 'phenotype' column: 0=control, 1=case)
    sample_metadata = bucket.schema("sample_metadata_schema").table("sample_metadata").select().read_all().to_pandas() # Get all metadata.
    merged_data = genotype_data.merge(sample_metadata[['sample_id', 'phenotype']], on='sample_id', how='left')
    # 3. Create a contingency table
    def genotype_to_counts(gt_str):
        if gt_str is None:
          return (0,0) # Missing
        gt = gt_str.replace('|', '/').split('/')
        return (gt.count('0'), gt.count('1'))


    merged_data['allele_counts'] = merged_data['gt'].apply(genotype_to_counts)

    merged_data[['ref_count', 'alt_count']] = merged_data['allele_counts'].apply(pd.Series)

    contingency_table = pd.crosstab(
        [merged_data['ref_count'],merged_data['alt_count']],
        merged_data['phenotype']
    )
    print(contingency_table)

    # 4. Perform the Chi-squared test
    chi2, p, dof, expected = stats.chi2_contingency(contingency_table)

    print(f"Chi-squared statistic: {chi2}")
    print(f"P-value: {p}")
    print(f"Degrees of freedom: {dof}")
    print(f"Expected frequencies: \n{expected}")

```

*   **Explanation:**
    *   This example demonstrates how to perform a Chi-squared test for association between a single variant and a binary phenotype (case/control status).
    *   It retrieves the genotype data for the variant and merges it with sample metadata (which is assumed to contain a `phenotype` column).
    *   It creates a contingency table using `pandas.crosstab()`.
    *   It uses `scipy.stats.chi2_contingency()` to perform the Chi-squared test and obtain the test statistic, p-value, degrees of freedom, and expected frequencies.
    *  **Important:** This example shows the basic steps.  For a real association study, you would need to:
        *   Loop through multiple variants.
        *   Correct for multiple testing (e.g., using Bonferroni correction or FDR).
        *   Consider population stratification and other confounding factors.

Okay, continuing from where we left off, here's the section on calculating Linkage Disequilibrium (LD), and the subsequent sections to complete Chapter 5:

**5.3.2. Calculating Linkage Disequilibrium (LD)**

LD calculation typically involves comparing the observed frequency of haplotypes (combinations of alleles at different loci) to the expected frequency under the assumption of independence. Here's a Python example using `scipy.stats` and data retrieved from VAST DB:

```python
import scipy.stats as stats
import pandas as pd
import numpy as np  # Import numpy

with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)
    genotypes_table = bucket.schema("genotypes_schema").table("genotypes")

    # Example: Calculate LD between two variants on chromosome 22
    pos1 = 17000050
    pos2 = 17000100

    # 1. Retrieve genotypes for both variants
    reader1 = genotypes_table.select(
        columns=['sample_id', 'gt'],
        predicate=(
            (_.chrom == "chr22") & (_.pos == pos1)
        )
    )
    genotypes1 = pa.Table.from_batches([b for b in reader1]).to_pandas()

    reader2 = genotypes_table.select(
        columns=['sample_id', 'gt'],
        predicate=(
            (_.chrom == "chr22") & (_.pos == pos2)
        )
    )
    genotypes2 = pa.Table.from_batches([b for b in reader2]).to_pandas()

    # 2. Merge the genotypes based on sample_id
    merged_genotypes = genotypes1.merge(genotypes2, on='sample_id', suffixes=('_1', '_2'))

    # 3. Convert genotypes to allele counts (0 and 1, assuming biallelic)
    def alleles(gt_str):
        if gt_str is None:
            return (None, None)  # Missing genotype
        return tuple(gt_str.replace('|', '/').split('/'))

    merged_genotypes['alleles_1'] = merged_genotypes['gt_1'].apply(alleles)
    merged_genotypes['alleles_2'] = merged_genotypes['gt_2'].apply(alleles)

    # 4. Create a contingency table of haplotype counts
    #    (This is a simplified example; you'd need to handle phasing and missing data more robustly)
    haplotype_counts = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for _, row in merged_genotypes.iterrows():
        a1_1, a2_1 = row['alleles_1']
        a1_2, a2_2 = row['alleles_2']

        #Handle missing data
        if a1_1 is None or a1_2 is None or a2_1 is None or a2_2 is None:
            continue
        
        # Count assuming haplotype (a1_1, a1_2)
        haplotype_counts[(int(a1_1), int(a1_2))] += 1
        # Count assuming haplotype (a2_1, a2_2)
        haplotype_counts[(int(a2_1), int(a2_2))] += 1


    # 5. Calculate D, D', and r^2
    n = sum(haplotype_counts.values())  # Total number of haplotypes
    if n == 0:
        print("No data to calculate LD.")
    else:
        p_A = (haplotype_counts[(0, 0)] + haplotype_counts[(0, 1)]) / n  # Frequency of allele 0 at locus 1
        p_B = (haplotype_counts[(0, 0)] + haplotype_counts[(1, 0)]) / n  # Frequency of allele 0 at locus 2
        p_AB = haplotype_counts[(0, 0)] / n  # Frequency of haplotype (0, 0)

        D = p_AB - (p_A * p_B)

        if D >= 0:
            D_max = min(p_A * (1 - p_B), (1 - p_A) * p_B)
        else:
            D_max = min(p_A * p_B, (1 - p_A) * (1 - p_B))

        D_prime = D / D_max if D_max != 0 else 0

        r_squared = D**2 / (p_A * (1 - p_A) * p_B * (1 - p_B)) if (p_A * (1 - p_A) * p_B * (1 - p_B)) !=0 else 0

        print(f"D: {D}")
        print(f"D': {D_prime}")
        print(f"r^2: {r_squared}")

```

*   **Explanation:**
    *   This example retrieves genotypes for two variants and merges them based on `sample_id`.
    *   It converts genotypes to allele counts.
    *   It constructs a contingency table of haplotype counts.
    *   It calculates D, D', and r^2, common measures of LD.
    *   **Important:** This is a simplified LD calculation.  A production-ready implementation would need to:
        *   Handle missing genotypes appropriately.
        *   Account for phasing (the example assumes unphased data).
        *   Potentially use more efficient algorithms for calculating LD across many pairs of variants.  Libraries like `scikit-allel` provide optimized functions for this.

**5.3.3. Integrating with External Statistical Packages (e.g., R, Python)**

The examples in 5.3.1 and 5.3.2 already demonstrate integration with Python's `scipy.stats` library. You can similarly integrate with other statistical packages:

*   **R:** Use the `reticulate` package in Python to call R functions and exchange data between Python and R.  You can retrieve data from VAST DB using the Python SDK, pass it to R for analysis, and then potentially store results back in VAST DB.
*  **Pandas/Numpy/Scipy**: As shown, you can readily use these libraries after retrieving data.

**5.4. Working with Snapshots (clones): consistent data views, rollbacks.**

VAST DB supports snapshots, which provide a consistent view of the data at a specific point in time.  This is extremely useful for:

*   **Reproducibility:**  Ensure that analyses are performed on a consistent dataset, even if the underlying data is being updated.
*   **Experimentation:** Create a snapshot, perform experimental analyses or data modifications, and then easily roll back to the original state if needed.
*   **Data Backup:**  Snapshots can serve as a form of data backup.

```python
with session.transaction() as tx:
    bucket = tx.bucket(DATABASE_NAME)

    # Create a snapshot of the entire bucket
    try:
        bucket.create_snapshot("my_snapshot")
        print("Snapshot 'my_snapshot' created.")
    except vastdb.errors.SnapshotExists:
        print("Snapshot already exists.")

    # List existing snapshots.
    print(bucket.snapshots())

with session.transaction() as tx:
     bucket = tx.bucket(DATABASE_NAME)
     # Access data from a snapshot (using a read-only transaction)
     #  You can use snapshot.bucket() to retrieve the tables *as of the snapshot*.

     try:
         snapshot= bucket.snapshot("my_snapshot")
         if snapshot:

            #Access the snapshot data
            variants_table_snap = snapshot.bucket().schema("variants_schema").table("variants")
            reader = variants_table_snap.select() # Select all data
            snapshot_data = pa.Table.from_batches([b for b in reader])
            print(snapshot_data.to_pandas().head())

            #Try to write to a table.
            try:
                variants_table_snap.insert([["FAIL", 1, "FAIL", "FAIL","FAIL"]]) #THIS WILL CAUSE AN EXCEPTION
            except Exception as e:
                print(f"Exception raised as expected: {e}")

            #Clean up the snapshot.
            snapshot.drop()
     except vastdb.errors.MissingSnapshot as e:
        print(f"Snapshot missing {e}")

```

*   **Explanation:**
    *   `bucket.create_snapshot("my_snapshot")`:  Creates a snapshot named "my_snapshot" of the entire bucket.
    *    `bucket.snapshots()`: List snapshots.
    *    `bucket.snapshot("my_snapshot")`: Get a handle to a named snapshot.
    *   `snapshot.bucket()`: Within a read-only transaction (note there is NO `with session.transaction() as tx` wrapping this block), you access the tables *as they existed at the time of the snapshot*. Any attempts to modify data within this context will result in an error. This ensures data consistency.
     *    `snapshot.drop()`: Delete a named snapshot.
* **Important considerations:**
    * You can not create a table within a shapshot context.
    * Be aware of the snapshot limit (1019) in VAST DB.
    * Snapshots consume storage space. Manage snapshots carefully, deleting them when they are no longer needed.
