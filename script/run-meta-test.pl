#! /usr/bin/perl

$list_file = $ARGV[0];
$type = $ARGV[1];
$prefix = $ARGV[2];

open(LIST, "<", $list_file);

for ($i = 0; $i < 5; ++$i)
{
	$line = <LIST>;
	@items = split(/\s+/, $line);
	$ref1 = $items[6];
	$ref2 = $items[7];
	$ref3 = $items[8];

	$ref = $prefix . "/ref.fa";
	$read = $prefix . "/read.fa";

	$read1 = $prefix . "/read1.fa";
	$read2 = $prefix . "/read2.fa";
	$read3 = $prefix . "/read3.fa";

	$cmd = "bin/sim_reads --paired --depth 10 $ref1 $read1";
	print $cmd, "\n";
	system($cmd);

	$cmd = "bin/sim_reads --paired --depth 100 $ref2 $read2";
	print $cmd, "\n";
	system($cmd);

	$cmd = "bin/sim_reads --paired --depth 1000 $ref3 $read3";
	print $cmd, "\n";
	system($cmd);

	$cmd = "cat $read1 $read2 $read3 > $read";
	print $cmd, "\n";
	system($cmd);

	$cmd = "cat $ref1 $ref2 $ref3 > $ref";
	print $cmd, "\n";
	system($cmd);

	$contig = "";
	$scaffold = "";
	if ($type eq "idba-ud")
	{
		$cmd = "time bin/idba_ud -r $read -o $prefix --pre_correction";
		print $cmd, "\n";
		system($cmd);

		$contig = $prefix . "/contig.fa";
		$scaffold = $prefix . "/scaffold.fa";
	}
	elsif ($type eq "soap")
	{
		$cmd = "cp soap-config $prefix";
		print $cmd, "\n";
		system($cmd);

		$cmd = "echo 'p=$read' >> $prefix/soap-config";
		print $cmd, "\n";
		system($cmd);

		$cmd = "time SOAPdenovo-127mer all -s $prefix/soap-config -K 45 -o $prefix/soap ";
		print $cmd, "\n";
		system($cmd);

		$contig = $prefix . "/soap.contig";
		$scaffold = $prefix . "/soap.scafSeq";
	}
	elsif ($type eq "velvet")
	{
		$cmd = "time script/velvet-pe $prefix/velvet 45 $read";
		print $cmd, "\n";
		system($cmd);

		$cmd = "bin/split_scaffold $prefix/velvet/contigs.fa";
		print $cmd, "\n";
		system($cmd);

		$contig = "$prefix/velvet/contigs.fa.split";
		$scaffold = "$prefix/velvet/contigs.fa";
	}
	elsif ($type eq "meta-idba")
	{
		$cmd = "time /home/pengyu/programming/release/idba-0.19/bin/metaidba -r $read -o $prefix/meta --mink 20 --maxk 100";
		print $cmd, "\n";
		system($cmd);

		$contig = $prefix . "/meta.consensus";
		$scaffold = $prefix . "/meta.consensus";
	}


	$cmd = "script/validate_blat_parallel $ref1 $contig 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref1 $scaffold 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref2 $contig 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref2 $scaffold 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref3 $contig 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref3 $scaffold 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref $contig 0.95 100";
	print $cmd, "\n";
	system($cmd);

	$cmd = "script/validate_blat_parallel $ref $scaffold 0.95 100";
	print $cmd, "\n";
	system($cmd);
}

