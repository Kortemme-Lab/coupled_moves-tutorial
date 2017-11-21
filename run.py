#!/usr/bin/python

import socket
import sys
import os
import subprocess

use_multiprocessing = False
if use_multiprocessing:
    import multiprocessing

coupled_moves_path = os.path.expanduser("~/rosetta/working_branches/ddg_backrub/source/bin/coupled_moves.linuxclangrelease")
nstruct = 1 # Would be 20 in normal usage

print "Python:", sys.version
print "Host:", socket.gethostname()

f = open('benchmark2.txt')
systems = []
for line in f:
    name = line.split()[0]
    extra = ''
    if len(line.split()) > 1:
        extra = line.split()[1]
    systems.append( [name, extra] )
f.close()

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

if use_multiprocessing:
    pool = multiprocessing.Pool()
for nstruct_i in range(1, nstruct + 1 ):
    for name, extra in systems:
        if use_multiprocessing:
            pool.apply_async( run_coupled_moves, args = (name, extra, nstruct_i) )
        else:
            run_coupled_moves( name, extra, nstruct_i )
if use_multiprocessing:
    pool.close()
    pool.join()
