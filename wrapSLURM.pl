use warnings;
use strict;
use Getopt::Long;
use Cwd; 

my ($commands, $nodes, $ppn, $mem, $walltime, $modules, $queue, $name, $tmpMem);

GetOptions( "commands=s" => \$commands,
		   "nodes=i" => \$nodes,
		   "ppn=i" => \$ppn,
		   "mem=s" => \$mem,
		   "tmpMem=s" => \$tmpMem,
		   "walltime=s" => \$walltime,
		   "modules=s" => \$modules,
		   "queue=s" => \$queue,
		   "jobname=s" => \$name);

$nodes = defined($nodes) ? $nodes : 1;
$ppn = defined($ppn) ? $ppn : 1;
$mem = defined($mem) ? $mem : "10gb";
$walltime = defined($walltime) ? $walltime : "6:00:00";
$queue = defined($queue) ? $queue : "ccgg";
$name = defined($name) ? $name : "wrappedJob";
$tmpMem = defined($tmpMem) ? $tmpMem : "10gb";

my $loadMods;
if (defined($modules)){
	open(MOD, "<", $modules) or die "can't open modules to load";
	while (<MOD>){
		my $line = $_;
		chomp($line);
		$loadMods = "$line ;";
	}
}
else{
	$loadMods = "";
}

my $sbatch = "sbatch -e ${name}.%j.e -o ${name}.%j.o -p $queue --tmp=$tmpMem --time=$walltime --ntasks-per-node=$ppn --mem=$mem --nodes=$nodes --job-name=$name";

open(COM, "<", $commands);
while(<COM>){
my $line = $_;
chomp($line);
print "$sbatch --wrap=\" $loadMods $line\"\n"; 
}