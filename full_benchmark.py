#!/usr/bin/env python2
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=48:00:00
#$ -t 1-250
#$ -j y
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -l netapp=1G

import sys
import socket
import os

print "Python:", sys.version
print "Host:", socket.gethostname()

sge_task_id = long( os.environ["SGE_TASK_ID"] )

coupled_moves_path = os.path.expanduser("~/bin/rosetta_static/coupled_moves.static.linuxgccrelease")
nstruct = 50
total_jobs = 250

def run_coupled_moves( name, extra, nstruct_i ):
    pdb = name.split('_')[0]
    lig = name.split('_')[1]
    pdb_file = os.path.join( name, pdb + "_with_" + lig + ".pdb" )
    params_file = os.path.join( name, lig + "_from_" + pdb + ".params" )
    res_file = os.path.join(name, pdb + '.resfile' )
    extra_params = ''
    if extra != '':
        extra_params = name + '/' + extra + "_from_" + pdb + ".params"

    output_directory = os.path.join( 'output', os.path.join( name, '%02d' % nstruct_i ) )
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    coupled_moves_args = [
        os.path.abspath(coupled_moves_path),
        "-s %s" % os.path.abspath(pdb_file),
        "-resfile %s" % os.path.abspath(res_file),
        "-extra_res_fa %s" % os.path.abspath(params_file),
        "-mute protocols.backrub.BackrubMover",
        "-ex1",
        "-ex2",
        "-extrachi_cutoff 0",
        "-nstruct %d" % 1,
        "-coupled_moves::mc_kt 0.6",
        "-coupled_moves::initial_repack false",
        "-coupled_moves::ligand_mode true",
        "-coupled_moves::fix_backbone false",
        "-coupled_moves::bias_sampling true",
        "-coupled_moves::boltzmann_kt 0.6",
        "-coupled_moves::bump_check true",
    ]

    if extra_params != '':
        coupled_moves_args.append("-extra_res_fa %s" % os.path.abspath(extra_params))

    log_path = os.path.join(output_directory, name + '.log')

    print 'Running Rosetta with args:'
    print ' '.join(coupled_moves_args)
    print 'Output logged to:', os.path.abspath(log_path)
    print

    outfile = open(log_path, 'w')
    process = subprocess.Popen(coupled_moves_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = output_directory)
    returncode = process.wait()
    outfile.close()

systems = []
with open('full_benchmark.txt', 'r') as f:
    for line in f:
        if len(line.strip()) > 0:
            name = line.split()[0]
            extra = ''
            if len(line.split()) > 1:
                extra = line.split()[1]
            systems.append( [name, extra] )

job_args = []
for nstruct_i in range(1, nstruct + 1 ):
    for name, extra in systems:
        job_args.append( (name, extra, nstruct_i) )
print nstruct, len(systems), len(job_args), total_jobs
assert( len(job_args) == total_jobs )

run_coupled_moves( job_args[sge_task_id-1] )
