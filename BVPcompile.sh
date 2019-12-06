#!/bin/bash
echo off
gcc ./Code/BVP.c ./Code/BVP_control.c ./Code/BVP_IC.c ./Code/BVP_sys.c ./Code/BVP_gridsize.c ./Code/BVP_out.c ./Code/BVP_physics.c ./Code/BVP_LSsolver.c ./Code/BVP_LUdcmp.c ./Code/BVP_BICO.c ./Code/BVP_GMRES.c ./Funcs/FEqEETTout.c ./Funcs/FEqEERRout.c ./Funcs/FEqEEPSIout.c ./Funcs/FEqISCOLRout.c -lm -o ./Execs/BVP.exe
echo BVP compiled
