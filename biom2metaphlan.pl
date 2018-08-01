#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Math::Round;
use Cwd 'abs_path';
use Data::Dumper;

my $VERSION = '1.0' ;
my $lastmodif = '2017-10-16' ;

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
) ;


&main;


sub main {
	my $self={};
	bless $self;
	$self->{_taxonIDs} = ["k__","p__","c__","o__","f__","g__", "s__"];
	_set_options($self);
	# $self->{_selected_category} = ['J3_actrade','J3_surgras','J3_rien'];
	$self->_read_map_file();
	$self->_read_biom_tab_file();
	my @taxo=();
	print 'sample_id' . "\t" . join("\t",keys(%{$self->{_map_hash}})) . "\n";
	foreach my $sample (keys(%{$self->{_map_hash}})){
		if(! defined( $self->{_map_hash}->{$sample} ) ){
			print STDERR $sample . "\n";
		}
		print "\t" . $self->{_map_hash}->{$sample}->{$self->{_category}};
	}
	print "\n";
	$self->_print($self->{_tree},\@taxo,-1);
}


sub _print {

	my ($self,$tree,$taxo,$depth)=@_;
	my %keysList;
	$depth++;
	# print Dumper $tree;
	foreach my $k (keys %{$tree}){
		if(!defined($self->{_sample_hash}->{$k})){
			$keysList{$k} = 1;
		}
	}
	foreach my $c (keys %keysList){
		my $newHash;
		if($self->{_is_taxonIDs} == 0){
			$taxo->[$depth] = $c;
		}
		else{
			$taxo->[$depth] = $self->{_taxonIDs}->[$depth] . $c;
		}
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
	open(CSV,$self->{_map_file}) || die('Cannot open file ' . $self->{_map_file});
	my $hash;
	my $header_line = <CSV>;
	chomp $header_line;
	my @headers = split("\t",$header_line);
	while(<CSV>){
		chomp;
		my @line = split("\t",$_);
		for(my $i=1;$i<=$#headers;$i++){
			# if(defined($self->{_selected_category})){
				if($headers[$i] eq $self->{_category}){
					# foreach my $s_cat (@{$self->{_selected_category}}){
						# if($line[$i] eq $s_cat){
							$hash->{$line[0]}->{$headers[$i]} = $line[$i];
						# }
					# }
				}
			# }
			else{
				$hash->{$line[0]}->{$headers[$i]} = $line[$i];
			}
		}
	}
	$self->{_map_hash} = $hash;
}


sub _read_biom_tab_file {
	my ($self)=@_;
	open(CSV,$self->{_input_file}) || die('Cannot open file ' . $self->{_input_file});
	my $tree ={};
	while(<CSV>){
		chomp;
		if(/^#\s/){
			next;
		}
		if(/^#OTU ID/){
			my @headers = split(/\t/,$_);
			for(my $i=1;$i<$#headers;$i++){
				if(defined($self->{_map_hash}->{$headers[$i]})){
					push(@{$self->{_sample_list}},$headers[$i]);
					$self->{_sample_hash}->{$headers[$i]}=$i;
				}
			}
		}
		else{
			my $parent = $tree;
			my @data_line = split(/\t/,$_);
			my @taxo = split(';',$data_line[$#data_line]);
			if(!defined($self->{_is_taxonIDs})){
				if($taxo[0] =~ /^k__/){
					$self->{_is_taxonIDs} = 1;
				}
				else{
					$self->{_is_taxonIDs} = 0;
				}
			}
			@taxo = $self->_treat_missing_rank(@taxo);
			foreach my $s (keys(%{$self->{_sample_hash}})){
				if(defined($self->{_map_hash}->{$s})){
					$self->{_total_per_sample}->{ $s } += $data_line[$self->{_sample_hash}->{$s}];
				}
			}
			for(my $j=0;$j<=$#taxo;$j++){
				$taxo[$j] =~ s/^\s//;
				foreach my $s (keys(%{$self->{_sample_hash}})){
					if(defined($self->{_map_hash}->{$s})){
						$parent->{$taxo[$j]}->{ $s } += $data_line[$self->{_sample_hash}->{$s}];
					}
				}
				$parent = $parent->{$taxo[$j]};
			}
		}
	}
	close CSV;
	$self->{_tree} = $tree;
}

sub _treat_missing_rank {
	my ($self,@taxo)=@_;
	my @new;
	my $last=0;
	for(my $i=0;$i<=$#taxo;$i++){
		my $tmp;
		if($self->{_is_taxonIDs} == 1){
			$tmp = $taxo[$i] ;
			$tmp =~ s/$self->{_taxonIDs}->[$i]//;
		}
		if($tmp eq ''){
			$last = $i;
			last;
		}
	}
	my $last_rank = $taxo[$last-1];
	$last_rank =~ s/$self->{_taxonIDs}->[$last-1]//;
	for(my $i=0;$i<=$#taxo;$i++){
		my $tmp;
		$tmp = $taxo[$i] ;
		if($self->{_is_taxonIDs} == 1){
			$tmp =~ s/$self->{_taxonIDs}->[$i]//;
		}
		if($i < $last){
			push(@new,$taxo[$i]);
		}
		else{
			push(@new,$self->{_taxonIDs}->[$i] . $last_rank);
		}
	}
	return @new;
}


sub _set_options {
	my ($self)=@_;
	if($input_file ne ''){
		$self->{_input_file} = abs_path($input_file);
	}
	else{
		print STDERR 'You must provide at least one tsv file.' . "\n";
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
	if($help){
		&help
	}
}



sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    This script generates Metaphlan compatible files from Biom TSV format.
#

USAGE: perl $prog -i otu-table-w-taxo.tsv -m mapping-file.tsv

				### OPTIONS ###
				-i|input        <OTU TABLE> OTU table in Biom TSV format, with taxonomy.
				-m|map          <QIIME MAP FILE> TSV file.
				-o|output       <FILE>  Output file name.
				-c              Map category to print.
				-d              Taxonomy depth. [int]
				-freq           Compute per sample frequencies.
				-help|h         Print this help and exit.
EOF
exit(1) ;
}
