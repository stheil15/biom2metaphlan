# biom2metaphlan


## Perl libraries

- Getopt::Long
- File::Basename
- Math::Round
- Cwd 'abs_path'

## Usage

```bash
biom2metaphlan.pl -i otu-table.tsv -m map.txt -c category -d [int]
```

## Help

```bash
### biom2metaphlan.pl 1.0 ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    1.0
# LAST MODIF: 2017-10-16
# PURPOSE:    This script is used to parse csv file containing tax_id field and creates Krona charts.
#

USAGE: perl biom2metaphlan.pl -i otu-table-w-taxo.tsv -m mapping-file.tsv

				### OPTIONS ###
				-i|input        <BLAST CSV>=<GROUP FILE>  Blast CSV extended file and CSV group file corresponding to blast (optional)
				-m|map          <QIIME MAP FILE> TSV file.
				-o|output       <DIRECTORY>  Output directory.
				-c              Map category to print.
				-d              Depth.
				-freq           Compute per sample frequencies.
				-help|h         Print this help and exit.
```
