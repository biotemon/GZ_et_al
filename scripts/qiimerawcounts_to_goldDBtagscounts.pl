#!/usr/bin/perl
#This script uses as input a table of qiime taxonomy tags and instert them into our gold taxonomic DB.
#Those tags are only added if they are not present in the db.
#Additionally, a table is generated as output with all the taxonomy columns of the DB along
#with its corresponding counts. 
#The script does two sets of steps A) check and update the db, and
#B) prepare and print the taxonomyXcounts table.
#Requirements: - Gold DB
#              - Longformat taxonomy tag x counts file

use DBI;
use strict;
use warnings;
use Term::ProgressBar;
use 5.010;
use List::MoreUtils qw(uniq);
use Data::Dumper qw(Dumper);

#SUBROUTINE TO CHECK SOMETHING IS IN AN ARRAY

sub in(&@){
  local $_;
  my $code = shift;
  for( @_ ){ # sets $_
    if( $code->() ){
      return 1;
    }
  }
  return 0;
}

#--OPEN FILES--

my ($driver, $database, $dsn, $userid, $password, $dbh);

#GZ_db
use DBI;
$driver = "SQLite";
#Next line should look like this.
#$database = '/Users/tito-admin/Tito/JOYELABACKUP/github/db/SK_metaT_v2.db';
$database = 'SETDABASEHERE';
$dsn = "DBI:$driver:dbname=$database"; #No spaces here!!

$userid = "";
$password = "";
$dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
or die $DBI::errstr;

my(@reca_blast, %counts_of_id, $val, @fields, @fields1, @fields1A, @fields3, @fields4, $id);
my ($i, $cmd, $outcmd, $progress_bar, @fields2, $stmt, $species, @input_file);

#Open list file
open INPUT, $ARGV[0];
@input_file = <INPUT>;
close INPUT;

my (@row, @row1, @add_manually_terms);

print "Determining taxonomic query terms\n";
$progress_bar = Term::ProgressBar->new(scalar(@input_file));
$i = 0;

my($cmd3, $cmd4, $outcmd4, $rsuperkingdom, $rkingdom, $rclass, $rorder);
my($rphylum, $rfamily, $rgenus, $rspecies, $itisid, @final_set, %taxid_of, @set);
my(@unique_query_terms);
@unique_query_terms = ();

#Check for unique terms
for (my $i=0; $i < scalar(@input_file); $i++)
{
	chomp ($input_file[$i]);
	@fields = split "\t", $input_file[$i];
	$fields[0] =~ s/\(//g;
 	$fields[0] =~ s/\)//g;
 	$fields[0] =~ s/\[//g;
 	$fields[0] =~ s/\]//g;
	if(!(in { $fields[0] eq $_ } @unique_query_terms )){
		push @unique_query_terms, $fields[0];
	 }
	$progress_bar->update($i);
}

print "Updating taxonomy database and getting tax_ids\n";
$progress_bar = Term::ProgressBar->new(scalar(@unique_query_terms));

foreach my $id (@unique_query_terms)
{
	$id =~ s/;__/;/g;
	$rkingdom = '';
	$rphylum = '';
	$rclass = '';
	$rorder = '';
	$rfamily = '';
	$rgenus = '';
	$rspecies = '';

	@fields1 = split(";", $id);

		if($fields1[0]){
		$rkingdom = $fields1[0]; 
		$rkingdom  =~ s/k_+//g;
		}

		if($fields1[1]){
		$rphylum = $fields1[1]; 
		$rphylum  =~ s/p_+//g;
		}

		if($fields1[2]){
		$rclass = $fields1[2]; 
		$rclass  =~ s/c_+//g;
		}

		if($fields1[3]){
		$rorder = $fields1[3]; 
		$rorder =~ s/o_+//g;
		}

		if($fields1[4]){
		$rfamily = $fields1[4]; 
		$rfamily =~ s/f_+//g;
		}

		if($fields1[5]){
		$rgenus = $fields1[5]; 
		$rgenus =~ s/g?_+//g;
		}

		if($fields1[6]){
		$rspecies = $fields1[6]; 
		$rspecies =~ s/s?_+//g;
		}

		@set = ($rkingdom, $rphylum, $rclass, $rorder, $rfamily, $rgenus, $rspecies);
		for (my $z=0; $z < scalar(@set); $z++)
			{
			if( $set[$z] eq ''){
				$set[$z] = '_';
				}
			}

	#print $id."\n";
	#trying to find this features in the DB
	$outcmd = undef;
	$cmd = "sqlite3 ".$database." \"SELECT * FROM TAXONOMY WHERE KINGDOM = \'".$set[0]."\' AND PHYLUM =  \'".$set[1]."\' AND CLASS = \'".$set[2]."\' AND ORDER_TAX = \'".$set[3]."\' AND FAMILY = \'".$set[4]."\' AND GENUS = \'".$set[5]."\' AND SPECIES = \'".$set[6]."\';\"";
	$outcmd = `$cmd`;
	if($outcmd){
		@fields2 = split '\|', $outcmd;
	 	#We found taxonomy using the qiime fields

	 	$taxid_of{$id}=$fields2[0];
	 }else{
		#check if it is in no_rank
	 	$outcmd = undef;
		$cmd = "sqlite3 ".$database." \"SELECT * FROM TAXONOMY WHERE NO_RANK = \'".$id."\';\"";
		$outcmd = `$cmd`;
	 	if($outcmd){
	 		@fields3 = split '\|', $outcmd;
	 		#We found taxonomy using the qiime fields
	 		$taxid_of{$id}=$fields3[0];
	 	}else{
	 	#insert qiime taxonomy tags in DB
	 	$stmt = $dbh->prepare('INSERT INTO TAXONOMY VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
	 	$rspecies = $set[5]."=".$set[6];
	 	$rspecies =~ s/=_//;
		$stmt->execute(undef, "GreenGenes", "_", "_", $set[0], $set[1], $set[2], $set[3], $set[4], $set[5] , $set[6], '_', $id) or die "Couldn't execute statement: " . $stmt->errstr;
		
		#Request the taxid of the last inserted row			
		$outcmd = undef;
		$cmd = "sqlite3 ".$database." \"select seq from sqlite_sequence where name= \'".'TAXONOMY'."\';\"";
		$outcmd = `$cmd`;
			
	    @fields4 = split '\|', $outcmd;
	    $taxid_of{$id}=$fields4[0];	
		}           
	 }
$i++;
$progress_bar->update($i);
} #line 86 foreach my $id (@unique_query_terms)

#Processing final table
print "Processing taxonomyXcounts table\n";
$progress_bar = Term::ProgressBar->new(scalar(@input_file));

my($cm4, $sample_id, $the_counts);
$val = 'GENE_ID'.'$\'\t\''.'TAX_ID'.'$\'\t\''.'ASSEMBLY_ID'.'$\'\t\''.'READ_COUNTS'.'$\'\t\''.'SOURCE'.'$\'\t\''.'ITIS_NUMBER'.'$\'\t\''.'SUPERKINGDOM'.'$\'\t\''.'KINGDOM'.'$\'\t\''.'PHYLUM'.'$\'\t\''.'CLASS'.'$\'\t\''.'ORDER_TAX'.'$\'\t\''.'FAMILY'.'$\'\t\''.'GENUS'.'$\'\t\''.'SPECIES'.'$\'\t\''.'SUBSPECIES'.'$\'\t\''.'NO_RANK';
$cmd4 = 'echo '.$val.'  > '.$ARGV[0].'_taxonomyXcounts.txt';
system($cmd4);

for (my $i=0; $i < scalar(@input_file); $i++)
{
	chomp ($input_file[$i]);
	@fields = split "\t", $input_file[$i];
	$fields[0] =~ s/;__/;/g;
	$fields[0] =~ s/\(//g;
 	$fields[0] =~ s/\)//g;
 	$fields[0] =~ s/\[//g;
 	$fields[0] =~ s/\]//g;

 	$id = $fields[0];
 	$sample_id = $fields[1];
 	$the_counts = $fields[2];
 	if(defined($taxid_of{$id})){
			$stmt = $dbh->prepare('SELECT * FROM TAXONOMY WHERE TAX_ID = ?');
			$stmt->execute($taxid_of{$id}) or die $DBI::errstr;
			@row = $stmt->fetchrow_array();

			for (my $z=0; $z < scalar(@row); $z++)
			    		{
				    		if(defined $row[$z]){
				    			print '';
				    		}else{
				   				$row[$z] = '';
				   			}
				   			}

			$val = '\''.$id.'\t'.$taxid_of{$id}.'\t'.$sample_id.'\t'.$the_counts.'\t'.$row[1].'\t'.$row[2].'\t'.$row[3].'\t'.$row[4].'\t'.$row[5].'\t'.$row[6].'\t'.$row[7].'\t'.$row[8].'\t'.$row[9].'\t'.$row[10].'\t'.$row[11].'\t'.'_'.'\'';
			$cmd4 = 'echo '.$val.'  >> '.$ARGV[0].'_taxonomyXcounts.txt';
			system($cmd4);
	}else{

		#$val = $id.'$\'\t\''.'NO_TAX_ID'.'$\'\t\''.$sample_id.'$\'\t\''.$the_counts.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.''.'$\'\t\''.'';
		#$cmd4 = 'echo '.$val.'  >> '.$ARGV[0].'_taxonomyXcounts.txt';
		#system($cmd4);
			$val = '\''.$id.'\t'.'NO_TAX_ID'.'\t'.$sample_id.'\t'.$the_counts.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.''.'\t'.'_'.'\'';
			$cmd4 = 'echo '.$val.'  >> '.$ARGV[0].'_taxonomyXcounts.txt';
			system($cmd4);
	}
	$progress_bar->update($i);

 }
