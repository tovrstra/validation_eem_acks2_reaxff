#!/usr/bin/env bash

echo "EEM Pair of water molecules"
qcalc eem ffield_eem s66/2701_01WaterWater100.xyz
echo
echo "EEM Pair of water molecules, qtot=1"
qcalc eem ffield_eem s66/2701_01WaterWater100.xyz -q 1
echo
echo "EEM Pair of water molecules, ionize first water, qtot=0"
qcalc eem ffield_eem s66/2701_01WaterWater100.xyz --constrain s66/2701_01WaterWater100.constraints.txt
echo
echo "EEM Pair of water molecules, ionize first water, qtot=1"
qcalc eem ffield_eem s66/2701_01WaterWater100.xyz -q 1 --constrain s66/2701_01WaterWater100.constraints.txt
echo

echo "EEM Orthorhomic example"
qcalc eem ffield_eem cod_orthorhombic/2200807.xyz
echo
echo "EEM Triclinic example"
qcalc eem ffield_eem cod_triclinic/1508741.xyz
echo
echo '------------------'

echo "ACKS2 Pair of water molecules"
qcalc acks2 ffield_acks2 s66/2701_01WaterWater100.xyz
echo
echo "ACKS2 Pair of water molecules, qtot=1"
qcalc acks2 ffield_acks2 s66/2701_01WaterWater100.xyz -q 1
echo
echo "ACKS2 Pair of water molecules, ionize first water, qtot=0"
qcalc acks2 ffield_acks2 s66/2701_01WaterWater100.xyz --constrain s66/2701_01WaterWater100.constraints.txt
echo
echo "ACKS2 Pair of water molecules, ionize first water, qtot=1"
qcalc acks2 ffield_acks2 s66/2701_01WaterWater100.xyz -q 1 --constrain s66/2701_01WaterWater100.constraints.txt
echo

echo "ACKS2 Orthorhomic example"
qcalc acks2 ffield_acks2 cod_orthorhombic/2200807.xyz
echo
echo "ACKS2 Triclinic example"
qcalc acks2 ffield_acks2 cod_triclinic/1508741.xyz
