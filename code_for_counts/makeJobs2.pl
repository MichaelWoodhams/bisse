#!/usr/bin/env perl
# Like makeJobs1.pl, but code instead of matching old runs now
# implements the restrictions in text:
# If mu0 > lambda0, code 2
# If mu1 > lambda1, code 3
# Else code 1. The 432 'code 1's are the 'good' scenarios.
# Can grep out the code 2 and code 3's if we want to only run code 1s.
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
    $code = "1-" . ($i+1);
    my ($lambda0,$lambda1) = split (" ",$lambda);
    if ($mu0 > $lambda0) {$code = "2-". ($i+1);}
    if ($mu1 > $lambda1) {$code = "3-". ($i+1);}
    
    print("[ -e treeOutput/$code-tree.txt  ] || /usr/bin/R < randomStartsAmazon.R --slave --args treeOutput/$code-tree.txt nodeOutput/$code-node.txt $lambda $mu0 $mu1 $q01 $q10 $code\n");
    # It is not possible to have code 2 and code 3 conditions apply at same time.
}
