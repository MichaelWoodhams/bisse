#!/usr/bin/env perl
my $lambda, $i, $mu0, $mu1, $q01, $q10, $code;
my @lambdaTable = ("0.5 1.5","1 1","1.5 0.5","0.2 1.8","1.8 0.2");
my @muTable = (0.01,0.25,0.5,0.8);
my @qTable = (0.01,0.05,0.1);
for $i (0..719) {
    $lambda = $lambdaTable[int($i/144)];
    $q10 =  $qTable[int($i/48) % 3];
    $q01 =  $qTable[int($i/16) % 3];
    $mu1 = $muTable[int($i/4 ) % 4];
    $mu0 = $muTable[    $i     % 4];
    if ($i < 3*144) {
	$code = "1-" . ($i+1);
    } elsif ($i < 4*144) {
	$code = "2-" . ($i+1);
    } else {
	$code = "3-" . ($i+1);
    }
#    print("$code # $lambda # $mu0 # $mu1 # $q01 # $q10 #\n");
    print("[ -e treeOutput/$code-tree.txt  ] || /usr/bin/R < randomStartsAmazon.R --slave --args treeOutput/$code-tree.txt nodeOutput/$code-node.txt $lambda $mu0 $mu1 $q01 $q10 $code\n");

}
