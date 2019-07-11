# EChO
Enhanced Chromatin Occupancy (EChO)

Scripts associated with EChO fragment size profile analysis of CUT&RUN data

## Citing EChO

If you use EChO in your work, please cite the following manuscript: 

Meers MP, Janssens DH, Henikoff S. (2019). "Pioneer factor-nucleosome binding events during differentiation are motif-encoded". Molecular Cell, in press.

## Usage

    EChO_1.0.sh <region bed file> <fragment bed file> ["foci" | "matrix"] <output prefix>

## Input descriptions

**Region bed file:** Refers to a bed file of enriched peaks from you CUT&RUN data, called with a peak caller such as SEACR (https://github.com/FredHutch/SEACR). Must contain at least three columns: chromosome, start, and end.

**Fragment bed file:** A bed file containing the forward and reverse read coordinates for every CUT&RUN fragment in your dataset. The fragment bed file can be generated from paired-end BAM files as follows: 
 
  1) Process you BAM file using the bedtools bamtobed utility (https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) with the -bedpe flag, which generates a 7-column file in which the first 6 columns correspond to the coordinates for the two mate pairs in the read.

  2) Filter the output file to throw out any lines for which the forward read chromosome is different than the reverse read chromosome (i.e. where column 1 does not equal column 4).

  3) Generate a new bed file from the resulting output by selecting the chromosome (column 1), the 5' coordinate of the forward read (column 2), and the 3' coordinate of the reverse read (column 6).

  4) Filter this output to throw out any read pairs for which the fragment length is greater than 1000, which are likely to represent incorrect mapping positions.
The result is a three-column bed file that represents every fragment with read pairs mapped within 1kb of each other.

**"foci"/"matrix":** Using "foci" runs EChO in foci mode, which calls all local minima (foci) in the EChO profiles derived from every entry in your region bed file. Using "matrix" runs EChO in matrix mode, which generates a matrix of base pair-resolution EChO fragment size values spanning a 400bp window for every entry in your region bed file. NOTE: For matrix mode, the region bed file **must** contain a fourth column that denotes the single base-pair focus within the region that will be used to center the matrix.

**Output prefix:** A prefix for naming of output files.


## Output descriptions

### Foci mode:

    <output prefix>.EChO.bed

This bed file reports every single base pair EChO focus in the input dataset. Output fields are as follows: chromosome, start, end, fragment size mean, fragment size coefficient of variation (standard deviation/mean), span value used

### Matrix mode:

    <output prefix>.EChO.matrix
 
This text file has 401 columns and a number of rows equal to the number of region bed file entries for which matrix generation succeeded. Each column contains an EChO average fragment size value corresponding to every base from 200 bp upstream to 200 bp downstream of the focus used for centering (column 4 in the region bed file).


## Example code

    EChO_1.0.sh peaks.bed fragments.bed foci myouput

Runs EChO in "foci" mode

    EChO_1.0.sh peaks.bed fragments.bed matrix mymatrix

Runs EChO in "matrix" mode
