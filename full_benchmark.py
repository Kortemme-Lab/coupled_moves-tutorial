#!/usr/bin/env python2
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=48:00:00
#$ -t 1-250
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -l netapp=1G

from run import run_coupled_moves

import sys
import socket
import os

print "Python:", sys.version
print "Host:", socket.gethostname()

sge_task_id = long( os.environ["SGE_TASK_ID"] )

nstruct = 50
total_jobs = 250

systems = []
with open('full_benchmark.txt', 'r') as f:
    for line in f:
        name = line.split()[0]
        extra = ''
        if len(line.split()) > 1:
            extra = line.split()[1]
        systems.append( [name, extra] )

job_args = []
for nstruct_i in range(1, nstruct + 1 ):
    for name, extra in systems:
        job_args.append( (name, extra, nstruct_i) )
assert( len(job_args) == total_jobs )

run_coupled_moves( job_args[sge_task_id-1] )
