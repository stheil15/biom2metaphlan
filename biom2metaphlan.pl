#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Logger::Logger;
use Math::Round;
use Cwd 'abs_path';

my $VERSION = '1.0' ;
my $lastmodif = '2017-09-29' ;

my $input_file= '';
my $map_file='';
my $output = '';
my $help;
my $verbosity = 1;
my $data;
my $freq=0;
my $category = '';
my $depth=0;

GetOptions(
          "i|input:s"      => \$input_file,
					"m|map:s"        => \$map_file,
					"c=s"            => \$category,
					"d=i"            => \$depth,
					"freq"           => \$freq,
          "o|output=s"     => \$output,
          "h|help"         => \$help,
          "v|verbosity=i"  => \$verbosity,
) ;


if($verbosity > 1){
  Logger::Logger->changeMode( $verbosity );
}


&main;


sub main {
  my $self={};
  bless $self;
  _set_options($self);
	$self->_read_map_file();
	$self->_read_biom_tab_file();
	$self->{_taxonIDs} = ["k__","p__","c__","o__","f__","g__", "s__"];
	my @taxo=();
	print 'sample_id' . "\t" . join("\t",@{$self->{_sample_list}}) . "\n";
	print $self->{_category};
	foreach my $sample (@{$self->{_sample_list}}){
		print "\t" . $self->{_map_hash}->{$sample}->{$self->{_category}};
	}
	print "\n";
	$self->_print($self->{_tree},\@taxo,-1);
}


sub _print {
	my ($self,$tree,$taxo,$depth)=@_;
	my %keysList;
	$depth++;
  foreach my $k (keys %{$tree}){
		if(!defined($self->{_sample_hash}->{$k})){
			$keysList{$k} = 1;
		}
	}
	# if($depth < $self->{_max_depth}){
		foreach my $c (keys %keysList){
			my $newHash;

			$taxo->[$depth] = $self->{_taxonIDs}->[$depth] . $c;
			if($#{$taxo} > 1){
				for( my $i = $depth+1 ; $i <= $#{$taxo} ; $i ++){
					delete $taxo->[$i];
				}
			}

			foreach my $k (keys %{$tree->{$c}}){
				if(defined($self->{_sample_hash}->{$k})){
				}
				else{
					$newHash->{$k} = $tree->{$c}->{$k};
				}
			}
			$self->{_old_depth} = $depth;
			if($depth <= $self->{_max_depth}){
				print join('|',@{$taxo});
				for(my $i=0;$i<=$#{$self->{_sample_list}};$i++){
					if($self->{_freq}){
						print "\t" . nearest(0.00001, $tree->{$c}->{$self->{_sample_list}->[$i]} / $self->{_total_per_sample}->{$self->{_sample_list}->[$i]});
					}
					else{
						print "\t" . $tree->{$c}->{$self->{_sample_list}->[$i]};
					}
				}
				print "\n";
			}
			_print($self,$newHash,$taxo,$depth);
		}
}


sub _read_map_file {
	my ($self)=@_;
	open(CSV,$self->{_map_file}) || $logger->logdie('Cannot open file ' . $self->{_map_file});
	my $hash;
	my $header_line = <CSV>;
	chomp $header_line;
	my @headers = split("\t",$header_line);
	while(<CSV>){
		chomp;
		my @line = split("\t",$_);
		for(my $i=1;$i<=$#headers;$i++){
			$hash->{$line[0]}->{$headers[$i]} = $line[$i];
		}
	}
	$self->{_map_hash} = $hash;
}

sub _read_biom_tab_file {
	my ($self)=@_;
	open(CSV,$self->{_input_file}) || $logger->logdie('Cannot open file ' . $self->{_input_file});
	my $tree ={};
	while(<CSV>){
		chomp;
		if(/^#\s/){
			next;
		}
		if(/^#OTU ID/){
			my @headers = split(/\t/,$_);
			for(my $i=1;$i<$#headers;$i++){
				push(@{$self->{_sample_list}},$headers[$i]);
				$self->{_sample_hash}->{$headers[$i]}=1;
			}
		}
		else{
			my $parent = $tree;
			my @data_line = split(/\t/,$_);
			my @taxo = split(';',$data_line[$#data_line]);
			for(my $i=1;$i<$#data_line;$i++){
				$self->{_total_per_sample}->{ $self->{_sample_list}->[$i-1] } += $data_line[$i];
			}
			for(my $j=0;$j<=$#taxo;$j++){
				$taxo[$j] =~ s/^\s//;
				for(my $i=1;$i<$#data_line;$i++){
					$parent->{$taxo[$j]}->{$self->{_sample_list}->[$i-1]} += $data_line[$i];
				}
				$parent = $parent->{$taxo[$j]};
			}
		}
	}
	close CSV;
	$self->{_tree} = $tree;
}


sub _set_options {
  my ($self)=@_;

  if($input_file ne ''){
		$self->{_input_file} = abs_path($input_file);
  }
  else{
    $logger->error('You must provide at least one csv file.');
    &help;
  }
	if($map_file ne ''){
		$self->{_map_file} = $map_file;
	}
	if($output ne ''){
		$self->{_output_file} = $output;
	}
	else{
		$self->{_output_file} = 'biom2metaphlan.res';
	}
	if ($category ne ''){
		$self->{_category} = $category;
	}
	if($depth!=0){
		$self->{_max_depth} = $depth;
	}
	$self->{_freq} = $freq;
}



sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    This script is used to parse csv file containing tax_id field and creates Krona charts.
#

USAGE: perl $prog -i otu-table-w-taxo.tsv -m mapping-file.tsv

       ### OPTIONS ###
       -i|input        <BLAST CSV>=<GROUP FILE>  Blast CSV extended file and CSV group file corresponding to blast (optional)
       -m|map          <QIIME MAP FILE> TSV file.
       -o|output       <DIRECTORY>  Output directory.
			 -c              Map category to print.
			 -d              Depth.
			 -freq           Compute per sample frequencies.
       -help|h         Print this help and exit.
EOF
exit(1) ;
}
