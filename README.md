# biom2metaphlan

Parse a TSV file format to Metaphlan format to run Lefse and/or Hclust scripts.

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

USAGE: perl biom2metaphlan.pl -i otu-table-w-taxo.tsv -m mapping-file.tsv

				### OPTIONS ###
				-i|input        <OTU TABLE> OTU table in Biom TSV format, with taxonomy.
				-m|map          <QIIME MAP FILE> TSV file.
				-o|output       <DIRECTORY>  Output directory.
				-c              Map category to print.
				-d              Taxonomy depth. [int]
				-freq           Compute per sample frequencies.
				-help|h         Print this help and exit.
```
