#!/bin/sh
mpirun -n 4 ../../build/pe_exec square.dmg square.smb out
paraview out/out.pvtu
