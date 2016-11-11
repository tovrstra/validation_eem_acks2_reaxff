#!/usr/bin/env bash

echo "Pair of water molecules"
./compute_charges.py eem ffield s66/2701_01WaterWater100.xyz
echo
echo "Pair of water molecules, qtot=1"
./compute_charges.py eem ffield s66/2701_01WaterWater100.xyz -q 1
echo
echo "Pair of water molecules, ionize first water, qtot=0"
./compute_charges.py eem ffield s66/2701_01WaterWater100.xyz --constrain s66/2701_01WaterWater100.constraints.txt
echo
echo "Pair of water molecules, ionize first water, qtot=1"
./compute_charges.py eem ffield s66/2701_01WaterWater100.xyz -q 1 --constrain s66/2701_01WaterWater100.constraints.txt
echo

echo "Orthorhomic example"
./compute_charges.py eem ffield cod_orthorhombic/2200807.xyz
echo
echo "Triclinic example"
./compute_charges.py eem ffield cod_triclinic/1508741.xyz
