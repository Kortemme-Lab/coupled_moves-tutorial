#!/usr/bin/python

import socket
import sys
import os
import subprocess

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

for name, extra in systems:
    pdb = name.split('_')[0]
    lig = name.split('_')[1]
    dir = name
    pdb_file = dir+'/'+pdb+"_with_"+lig+".pdb"
    params_file = dir+'/'+lig+"_from_"+pdb+".params"
    res_file = dir+'/'+pdb+'.resfile'
    extra_params = ''
    if extra != '':
        extra_params = dir+'/'+extra+"_from_"+pdb+".params"


    coupled_moves_args = [
        coupled_moves_path,
        "-s %s" % pdb_file,
        "-resfile %s" % res_file,
        "-extra_res_fa %s" % params_file,
        "-mute protocols.backrub.BackrubMover",
        "-ex1",
        "-ex2",
        "-extrachi_cutoff 0",
        "-nstruct %d" % nstruct,
        "-coupled_moves::mc_kt 0.6",
        "-coupled_moves::initial_repack false",
        "-coupled_moves::ligand_mode true",
        "-coupled_moves::fix_backbone false",
        "-coupled_moves::bias_sampling true",
        "-coupled_moves::boltzmann_kt 0.6",
        "-coupled_moves::bump_check true",
    ]

    if extra_params != '':
        coupled_moves_args.append("-extra_res_fa %s" % extra_params)

    print 'Running Rosetta with args:'
    print ' '.join(coupled_moves_args)
    print 'Output logged to:', name + '.log'
    print

    outfile = open(name+'.log', 'w')
    process = subprocess.Popen(coupled_moves_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True)
    returncode = process.wait()
    outfile.close()
