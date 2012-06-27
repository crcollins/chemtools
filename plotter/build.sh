#!/bin/bash

WINPYDIR="c:/Python27"
WINPYTHON="wine $WINPYDIR/python.exe"
WINPWD=`winepath -w \`pwd\``
PYINSTALLER="/home/chris/Desktop/chem/pyinstaller"
PROGRAM=${PWD##*/}

$WINPYTHON $PYINSTALLER/pyinstaller.py --onefile --upx --tk $PROGRAM.py
python $PYINSTALLER/pyinstaller.py --onefile --upx --tk $PROGRAM.py

cp dist/$PROGRAM.exe ../bin/$PROGRAM/$PROGRAM.exe
cp dist/$PROGRAM ../bin/$PROGRAM/$PROGRAM
cp $PROGRAM.py ../bin/$PROGRAM/$PROGRAM.py

rm warn$PROGRAM.txt
rm logdict*.log
rm $PROGRAM.spec
rm tk*

rm -rf build
rm -rf dist
