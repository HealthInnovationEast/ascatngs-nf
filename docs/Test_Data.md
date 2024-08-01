# Test Data

Testing of this workflow can be completed using data provided by the authors of dockstore-cgpwgs.

## Reference files

Download and extract the required files:

```
cd $TESTDATA
mkdir GRCh37
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz | tar -C GRCh37 --strip-components 1 -zx
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/qcGenotype_GRCh37d5.tar.gz | tar -C GRCh37 --strip-components 1 -zx
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz | tar --no-anchor -C GRCh37 --strip-components 1 -zx SnpGcCorrections.tsv
```

## COLO-829 & COLO-829-BL

Low depth copies of this cell line are functional.  Workflow supports CRAM so minimisation of space is considered.

```
cd $TESTDATA

wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829.bam
samtools view --write-index -T GRCh37/genome.fa -Co GRCh37/COLO-829.cram -@ 7 --output-fmt-option normal COLO-829.bam
rm COLO-829.bam

wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829-BL.bam
samtools view --write-index -T GRCh37/genome.fa -Co GRCh37/COLO-829-BL.cram -@ 7 --output-fmt-option normal COLO-829-BL.bam
rm COLO-829-BL.bam
```
