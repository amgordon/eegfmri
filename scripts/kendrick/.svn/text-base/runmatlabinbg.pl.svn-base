#!/usr/bin/perl -w

# call like this:
#   /usr/bin/nohup runmatlabinbg.pl <flag> <command> <outputfile> > /dev/null &
# and then you can safely kill the ssh connection without disturbing the MATLAB process.
#
# we run MATLAB in the background and issue <command>.
# if <flag> is 0, we set the number of computational threads to 1 (-singleCompThread).
# if <flag> is 1, we allow multiple computational threads.
# the output and errors from <command> are written to <outputfile>.
# upon completion, <outputfile> is renamed "<outputfile>.done".
#
# note: edit this file with the specific MATLAB call that you want.
# the default MATLAB call that we use is:
#   /usr/local/matlab/r2009b/bin/matlab -nosplash -nodesktop -nodisplay -nojvm

use strict;
$| = 1;

if (scalar(@ARGV)!=3) {
  print "wrong number of arguments, so exiting\n";
  exit;
}

# delete existing .done file
`rm -f $ARGV[2].done`;

# put a record of the command in the file
`echo "$ARGV[1]" > $ARGV[2]`;

# record date to the file
`date >> $ARGV[2]`;

# run matlab and output to the file
my $cmd = "fprintf('now issuing command $ARGV[1] on host %s (PID %d).\\n\\n',gethostname,getpid); $ARGV[1]; fprintf('\\ncommand $ARGV[1] completed.\\n'); quit;";
if ($ARGV[0]==0) {
  `echo "$cmd" | /usr/local/matlab/r2009b/bin/matlab -nosplash -nodesktop -nodisplay -singleCompThread >> $ARGV[2] 2>&1`;   # -nojvm
} else {
  `echo "$cmd" | /usr/local/matlab/r2009b/bin/matlab -nosplash -nodesktop -nodisplay >> $ARGV[2] 2>&1`;
}

# record date to the file
`date >> $ARGV[2]`;

# rename to indicate completion
`mv $ARGV[2] $ARGV[2].done`;




#`/usr/local/matlab/r2009b/bin/matlab -singleCompThread -nosplash -nodesktop -nodisplay -nojvm -r "$cmd" >> $ARGV[1] 2>&1`;
# maxNumCompThreads(1); 
