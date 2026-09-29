# Minimal FASTQ Trimming Pipeline

A streamlined Nextflow pipeline for fastq trimming with AWS Batch support, extracted from NeonDisco.

## Directory Structure

```
.
├── main.nf                 # Main workflow file
├── nextflow.config         # Configuration file
├── modules/
│   └── trimming.nf        # Trimming process module
└── manifest.tsv           # Optional: sample manifest file
```

## Requirements

- Nextflow >= 23.04.0
- Docker (for local execution)
- AWS CLI configured (for AWS Batch)
- AWS Batch compute environment and job queue set up

## Quick Start

### Local Execution

Using input directory:
```bash
nextflow run main.nf \
    -profile local \
    --inputDir /path/to/fastq/files \
    --inputSource local \
    --outputDir ./results
```

Using manifest file:
```bash
nextflow run main.nf \
    -profile local \
    --manifestPath manifest.tsv \
    --inputSource local \
    --outputDir ./results
```

### AWS Batch Execution

**Important**: Before running on AWS Batch, edit `nextflow.config` to set:
- `process.queue` to your AWS Batch queue name
- `aws.region` to your AWS region

```bash
nextflow run main.nf \
    -profile awsbatch \
    --manifestPath s3://my-bucket/manifest.tsv \
    --inputSource s3 \
    --outputDir s3://my-bucket/results
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--inputDir` | null | Local directory containing paired FASTQ files |
| `--manifestPath` | null | Tab-delimited manifest file (TSV) |
| `--inputSource` | `local` | Input source: `local` or `s3` |
| `--outputDir` | `./outputs` | Output directory (supports S3 URIs) |

## Input Requirements

### Directory Input
- FASTQ files must follow naming pattern: `*_R1*.fastq.gz` and `*_R2*.fastq.gz`
- Files should be paired-end reads

### Manifest Input
Tab-delimited file with header:
```
sampleName	read1Path	read2Path
Sample001	/path/to/Sample001_R1.fastq.gz	/path/to/Sample001_R2.fastq.gz
```

## Output Structure

```
outputs/
├── trimmed_reads/
│   ├── Sample001_R1_trimmed.fastq.gz
│   ├── Sample001_R2_trimmed.fastq.gz
│   ├── Sample001_fastp.html
│   ├── Sample001_fastp.json
│   └── ...
└── reports/
    ├── execution_report.html
    ├── timeline.html
    ├── trace.txt
    └── flowchart.html
```

## Configuration Notes

### AWS Batch Setup

Before using AWS Batch, ensure you have:

1. **Compute Environment**: Created with appropriate instance types
2. **Job Queue**: Linked to your compute environment
3. **IAM Roles**: Proper permissions for Batch, S3, ECR, and CloudWatch
4. **S3 Buckets**: For input data and results

Update `nextflow.config`:
```groovy
awsbatch {
    process {
        executor = 'awsbatch'
        queue = 'your-batch-queue-name'  // Your queue name here
    }
    
    aws {
        region = 'us-east-1'  // Your region here
    }
}
```

### Trimming Parameters

Default fastp settings in `modules/trimming.nf`:
- Adapter detection: automatic for paired-end
- Quality threshold: Q20
- Minimum length: 50 bp
- Threads: 4 (configurable via `process.cpus`)

To modify trimming parameters, edit the `script` section in `modules/trimming.nf`.

## Resource Configuration

Adjust resources in `nextflow.config`:

```groovy
process {
    withName: 'TRIM_READS_FASTP' {
        cpus = 8          // Increase for faster processing
        memory = '16 GB'  // Increase for large files
        time = '4.h'
    }
}
```

## Troubleshooting

### "No such file" errors
- Verify FASTQ naming patterns match `*_R1*.fastq.gz` and `*_R2*.fastq.gz`
- Check file paths in manifest are absolute or relative to execution directory

### AWS Batch queue not found
- Update `process.queue` in `nextflow.config` with your actual queue name
- Verify queue exists: `aws batch describe-job-queues`

### Docker permission denied (local)
- Ensure Docker is running and user has permissions
- On Linux: add user to docker group or use sudo

### S3 access denied
- Check IAM role/user has S3 read/write permissions
- Verify AWS credentials are configured: `aws s3 ls`

## Differences from NeonDisco

This minimal pipeline extracts only:
- ✅ FASTQ trimming with fastp
- ✅ AWS Batch execution support
- ✅ Manifest and directory input modes

Removed from original NeonDisco:
- ❌ STAR alignment
- ❌ Fusion calling (Arriba, FusionCatcher, STAR-Fusion)
- ❌ HLA typing
- ❌ Neoantigen prediction
- ❌ All downstream analysis modules

## License

Inherits license from original NeonDisco repository.
