# VCF Exploration

## View Headers

```
bcftools view -h ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```
Outputs (truncated for brevity):

```
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20150218
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##source=1000GenomesPhase3Pipeline
##contig=<ID=1,assembly=b37,length=249250621>
##contig=<ID=2,assembly=b37,length=243199373>
...
##contig=<ID=20,assembly=b37,length=63025520>
##contig=<ID=21,assembly=b37,length=48129895>
##contig=<ID=22,assembly=b37,length=51304566>
##contig=<ID=GL000191.1,assembly=b37,length=106433>
##contig=<ID=GL000192.1,assembly=b37,length=547496>
...
##contig=<ID=GL000248.1,assembly=b37,length=39786>
##contig=<ID=GL000249.1,assembly=b37,length=38502>
##contig=<ID=MT,assembly=b37,length=16569>
##contig=<ID=NC_007605,assembly=b37,length=171823>
##contig=<ID=X,assembly=b37,length=155270560>
##contig=<ID=Y,assembly=b37,length=59373566>
##contig=<ID=hs37d5,assembly=b37,length=35477943>
##ALT=<ID=CNV,Description="Copy Number Polymorphism">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
##ALT=<ID=INS:ME:LINE1,Description="Insertion of LINE1 element">
##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">
##ALT=<ID=INS:MT,Description="Nuclear Mitochondrial Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CN0,Description="Copy number allele: 0 copies">
##ALT=<ID=CN1,Description="Copy number allele: 1 copy">
##ALT=<ID=CN2,Description="Copy number allele: 2 copies">
...
##ALT=<ID=CN123,Description="Copy number allele: 123 copies">
##ALT=<ID=CN124,Description="Copy number allele: 124 copies">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CS,Number=1,Type=String,Description="Source call set.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MC,Number=.,Type=String,Description="Merged calls.">
##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END<POLARITY; If there is only 5' OR 3' support for this call, will be NULL NULL for START and END">
##INFO=<ID=MEND,Number=1,Type=Integer,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Integer,Description="Estimated length of mitochondrial insert">
##INFO=<ID=MSTART,Number=1,Type=Integer,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length. It is only calculated for structural variation MEIs. For other types of SVs; one may calculate the SV length by INFO:END-START+1, or by finding the difference between lengthes of REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=TSD,Number=1,Type=String,Description="Precise Target Site Duplication for bases, if unknown, value will be NULL">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth; only low coverage data were counted towards the DP, exome data were not used">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele, IndelType:Type of Indel (REF, ALT and IndelType are only defined for indels)">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
##INFO=<ID=EX_TARGET,Number=0,Type=Flag,Description="indicates whether a variant is within the exon pull down target boundaries">
##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag,Description="indicates whether a site is multi-allelic">
##bcftools_viewVersion=1.21+htslib-1.21
##bcftools_viewCommand=view -h ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz; Date=Thu Mar 20 13:31:10 2025
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102	...
```

The last lines has 2512 items.  Therefore, the VCF file has 2512 columns of data.

```
liststr = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT HG00096	HG00097	HG00099	... "
l = liststr.split("\t")
print(len(l))
# 2512
```

## View some data

```
bcftools view ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | grep -E "^#CHROM|^[^#]" | head -n 2
```

Outputs:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099 ...
22	16050075	.	A	G	100	PASS	AC=1;AF=0.000199681;AN=5008;NS=2504;DP=8012;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0.001;AA=.|||;VT=SNP	GT	0|0	0|0	0|0 ...
```

## View specific columns

```
bcftools view -H -v snps ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | head -n 2 | cut -f 1-12
```

Outputs:

```
22	16050075	.	A	G	100	PASS	AC=1;AF=0.000199681;AN=5008;NS=2504;DP=8012;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0.001;AA=.|||;VT=SNP	GT	0|0	0|0	0|0
22	16050115	.	G	A	100	PASS	AC=32;AF=0.00638978;AN=5008;NS=2504;DP=11468;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0.0234;EUR_AF=0;SAS_AF=0;AA=.|||;VT=SNP	GT	0|0	0|0	0|0
```

## Description

Let's break down the first data record from the output of `bcftools view -H -v snps ... | head -n 2 | cut -f 1-12`:

```
22	16050075	.	A	G	100	PASS	AC=1;AF=0.000199681;AN=5008;NS=2504;DP=8012;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0.001;AA=.|||;VT=SNP	GT	0|0	0|0	0|0
```

**Column-by-Column Explanation:**

1.  **`22` (CHROM):**
    * This indicates the variant is located on chromosome 22.
2.  **`16050075` (POS):**
    * This is the position of the variant on chromosome 22, specifically at base position 16,050,075.
3.  **`.` (ID):**
    * This is an identifier for the variant. In this case, it's a dot (`.`), meaning no specific ID (like an rsID) is available.
4.  **`A` (REF):**
    * This is the reference allele. It means that the reference genome has an "A" (adenine) at this position.
5.  **`G` (ALT):**
    * This is the alternate allele. It means that at least one sample in the dataset has a "G" (guanine) at this position instead of an "A".
6.  **`100` (QUAL):**
    * This is the Phred-scaled quality score, indicating the confidence in the variant call. A higher score means higher confidence. 100 is a relatively high score.
7.  **`PASS` (FILTER):**
    * This indicates that the variant passed all filters.
8.  **`AC=1;AF=0.000199681;AN=5008;NS=2504;DP=8012;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0.001;AA=.|||;VT=SNP` (INFO):**
    * This is a semicolon-separated list of variant information:
        * `AC=1`: Allele count (1 alternate allele observed).
        * `AF=0.000199681`: Allele frequency (approximately 0.0002).
        * `AN=5008`: Total number of alleles in the dataset.
        * `NS=2504`: Number of samples with data.
        * `DP=8012`: Total read depth at this position.
        * `EAS_AF=0`, `AMR_AF=0`, `AFR_AF=0`, `EUR_AF=0`, `SAS_AF=0.001`: Allele frequencies in different populations (East Asian, Admixed American, African, European, South Asian).
        * `AA=.|||`: Ancestral allele information (unavailable).
        * `VT=SNP`: Variant type (single nucleotide polymorphism).
9.  **`GT` (FORMAT):**
    * This indicates that the genotype information is provided in the following columns.
10. **`0|0` (HG00096, HG00097, HG00099):**
    * These columns represent the genotypes for samples HG00096, HG00097, and HG00099 respectively.
    * `0|0` means the sample is homozygous for the reference allele.
        * The `0` refers to the reference allele.
        * The `|` indicates that the genotype is phased. Phasing means that the alleles are assigned to specific chromosomes.
        * If the `|` was a `/` it would mean the genotype is unphased.
    * **Genotype Interpretations:**
        * `0|0` or `0/0`: Homozygous reference (two copies of the reference allele).
        * `0|1` or `0/1`: Heterozygous (one copy of the reference allele, one copy of the alternate allele).
        * `1|0` or `1/0`: Heterozygous (same as above, but phased differently).
        * `1|1` or `1/1`: Homozygous alternate (two copies of the alternate allele).

**In Summary:**

The first record describes a SNP at position 16,050,075 on chromosome 22, where an "A" is replaced by a "G" in at least one sample. This variant is rare in the 1000 Genomes dataset. The first three samples, HG00096, HG00097, and HG00099, are all homozygous for the reference allele (A/A).

