#!/bin/bash

mkdir DEPENDENCIES 2>/dev/null
cd DEPENDENCIES

echo "******************************** BEGIN CUBA BUILD ********************************"
echo "******************************** BEGIN CUBA BUILD ********************************"
echo "******************************** BEGIN CUBA BUILD ********************************"
echo "******************************** BEGIN CUBA BUILD ********************************"
echo "******************************** BEGIN CUBA BUILD ********************************"
rm -r Cuba-3.2* 2> /dev/null
wget http://www.feynarts.de/cuba/Cuba-3.2.tar.gz
tar -xvzf Cuba-3.2.tar.gz
cd Cuba-3.2
sed -e "s/qmake/qmake-qt4/g" -i makefile.in
./configure  --prefix=/usr  --datadir=/usr/share/doc/ CFLAGS=" -fPIC -pipe -fstack-protector --param=ssp-buffer-size=4 -D_FORTIFY_SOURCE=2   -O3 " CPPFLAGS=" -pipe -fstack-protector --param=ssp-buffer-size=4 -D_FORTIFY_SOURCE=2  -fPIC -O3 "
mkdir ../Cuba-install
make install DESTDIR="../Cuba-install"
cd ..
rm -r Cuba-3.2
mv Cuba-install/usr Cuba-3.2
rm -r Cuba-install

echo ""
echo ""
echo ""
echo ""
if ls Cuba-3.2/include/cuba.h 2>/dev/null 1>/dev/null; then
    echo "Building Cuba was successful!"
    echo "CUBALIBS=../DEPENDENCIES/Cuba-3.2/lib/" >> ../install.env
    echo "CUBAINCLUDES=../DEPENDENCIES/Cuba-3.2/include/" >> ../install.env
else
    echo "Encountered a problem while building Cuba!"
    echo "This script is not smart enough fix it :("
fi
