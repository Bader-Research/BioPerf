# assumes that the environment variable THREADS has been set
# if not, you can do this with setenv THREADS <number>
# or export THREADS=<number>
# recommended is that <number> is one more than the number of CPUs on the
# machine

import os, sys, errno, time
from os.path import *
from stat import ST_SIZE

def stress(seqdir, smpExe, origExe, resultsFile):
    # Stress test smp and original version
    # and dump the results into a file
    #works on all .faa files in the seqdir directory
    try:
        # always open destructively
        rf = open(resultsFile, 'w')
        try:
            os.mkdir("smp_results")
            os.mkdir("orig_results")
        except:
            # directories already exist
            # it would be nice if they were empty but we dont care really
            pass
        absSmpExe = abspath(smpExe)
        absOrigExe = abspath(origExe)
        for seqfile in os.listdir(seqdir):
            name, ext = splitext(seqfile)
            if ext == '.faa' or ext == '.seq' or ext == '.fasta':
                # for now work only on files with '.faa' or '.seq' extension
                # take the starting time first for smp version
                absfile = seqdir+'/'+seqfile
                print absfile
                start = time.time()
                os.system("%s -infile=%s > /dev/null" %(absSmpExe,absfile))
                end = time.time()
                # now write the results
                rf.write('file: '+name+ext+'\n')
                rf.write('smp time: ' + str(end-start) + '\n')
                # move the resulting .dnd and .aln files to smp_results
                os.rename(seqdir+'/'+name+'.dnd', "smp_results/"+name+'.dnd')
                os.rename(seqdir+'/'+name+'.aln', "smp_results/"+name+'.aln')
                # now do the test on the original executable
                start = time.time()
                os.system("%s -infile=%s > /dev/null" %(absOrigExe,absfile))
                end = time.time()
                # now write the results
                rf.write('original time: ' +  str(end-start) + '\n')
                # move the resulting .dnd and .aln files to smp_results
                os.rename(seqdir+'/'+name+'.dnd', "orig_results/"+name+'.dnd')
                os.rename(seqdir+'/'+name+'.aln', "orig_results/"+name+'.aln')
                # and now run diff on the files
                os.system('diff smp_results/%s.dnd orig_results/%s.dnd > %s.dndiffs' %((name,)*3))
                os.system('diff smp_results/%s.aln orig_results/%s.aln > %s.alndiffs' %((name,)*3))
                if (os.stat(name+'.dndiffs')[ST_SIZE] == 0) and (os.stat(name+'.alndiffs')[ST_SIZE] == 0):
                    rf.write('Passed: OK\n')
                    os.unlink(name+'.dndiffs')
                    os.unlink(name+'.alndiffs')
                else:
                    rf.write('Failed\n')
        rf.close()
    except OSError:
        print sys.exc_type, sys.exc_value

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'python stress.py seqdir smpExe origExe resultFile <THREADS>'
        print 'results will be stored in resultFile'
        print 'clustalw_smp stress testing by Ognen Duzlevski'
        print '(C) 2001 Plant Biotechnology Institute, NRC'
        sys.exit(1)
    if len(sys.argv) == 6:
        # THREADS is given to us in command line
        os.putenv('THREADS', sys.argv[5])
    elif len(sys.argv) == 5:
        # check if THREADS is already set
        if not os.environ.has_key('THREADS'):
            print 'Set the THREADS environment variable.'
            print 'Use either setenv THREADS <number> or export THREADS=<number>'
            print 'Recommended number of threads is one more than CPUs'
            sys.exit(1)
    # do the stress testing
    stress(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
