#!/bin/bash
echo off
echo Watch me delete this shit...
rm -r "../../Data/BVPout_sols"
mkdir "../../Data/BVPout_sols"
mkdir "../../Data/BVPout_sols/Sol_Funcs"
mkdir "../../Data/BVPout_sols/Sol_Props"
echo Successfully deleted that shit.
rem BVP.exe 130 -5 101 12 1.0
BVP.exe 1 -5 101 12 1.0
FOR /L %%A IN (10,10,130) DO (
	BVP.exe %%A -5 101 12 1.0
)
echo Succesfully ran BVP.
