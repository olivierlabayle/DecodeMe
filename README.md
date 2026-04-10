# DecodeMe

This repository holds the instructions for reproducing the chromosome X DecodeME results.

## Notes / Questions

- call rate < 0.98 after merge --> If a variant has a bad call rate in DecodeMe (e.g. 5000 samples) it will still proceed. Should I enforce this for every batch?
- Was the imputation performed on the intersection of variants betweem DecodeMe and UKB ? Otherwise imputation might still be slightly biased.
- Is there a UKB control sample file? For now I have extracted one from the imputed chr1 fam file.
- To Dominique set parental id to 0 is not known ?

23       AX-94364057    T    -    0.0001168     8565 --> A1 and A2, total < 200 SNPs
23       AX-94364058    0    G            0     8564 --> A1 only, ~1200 SNPs in batch 1

## Prerequisites

- Install the [dx toolkit and upload agent](https://documentation.dnanexus.com/downloads) on your local machine.

- Install the [dxCompiler](https://github.com/dnanexus/dxCompiler?tab=readme-ov-file#setup).

- Login

```bash
dx login
```

I assume that both `dx` and `ua` are available in your path an that the `AUTH_TOKEN` and `RAP_PROJECT_ID` are set.

## Upload the DecodeME Data

The `DECODEME_DIRPATH` variables points to the folder where the genotype batches are stored.

```bash
DECODEME_DIRPATH=/gpfs/igmmfs01/datastore/DecodeME-DNA/chrX-analysis/chrX-Array-QC/for-imputation/first-test/
ua \
--auth-token ${AUTH_TOKEN} \
--project ${RAP_PROJECT_ID} \
--do-not-compress \
--progress \
--folder "/decodeme-chrX-genotyped" \
${DECODEME_DIRPATH}
```

## Merge with UK Biobank

First compile the workflow, `DX_COMPILER_PATH` points to the dxCompiler jar file.

```bash
java -jar $DX_COMPILER_PATH compile merge_with_ukb.wdl \
-f -project $RAP_PROJECT_ID \
-extras config/extras.json \
-reorg \
-folder /workflows/merge_with_ukb \
-inputs config/merge_with_ukb.json
```

Then run it:

```bash
dx run -y \
-f config/merge_with_ukb.dx.json \
--priority high \
--preserve-job-outputs \
--destination /chrX-merge-with-ukb/ \
/workflows/merge_with_ukb/merge_with_ukb
```

## Tips

To access a virtual machine and run code locally:

```bash
dx run \
--instance-type mem2_ssd1_v2_x8 \
-imax_session_length="10h" \
-y \
--ssh app-cloud_workstation
```