#!/bin/sh


SOURCE_DIR='STG_core/lib'
BUILD_DIR='build'

mkdir $BUILD_DIR
cp $SOURCE_DIR/*.c $BUILD_DIR
cp $SOURCE_DIR/*.h $BUILD_DIR
cp Makefile $BUILD_DIR

cd $BUILD_DIR

make -f Makefile clean
make -f Makefile

rm -f *.c
rm -f *.h

cd ..


