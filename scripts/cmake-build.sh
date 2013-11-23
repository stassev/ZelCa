#!/bin/bash

mkdir DEPENDENCIES 2>/dev/null
cd DEPENDENCIES

echo "******************************** BEGIN CMAKE BUILD ********************************"
echo "******************************** BEGIN CMAKE BUILD ********************************"
echo "******************************** BEGIN CMAKE BUILD ********************************"
echo "******************************** BEGIN CMAKE BUILD ********************************"
echo "******************************** BEGIN CMAKE BUILD ********************************"
rm -r cmake-2.8.12.1* 2> /dev/null
wget http://www.cmake.org/files/v2.8/cmake-2.8.12.1.tar.gz
tar -xzvf cmake-2.8.12.1.tar.gz
cd cmake-2.8.12.1
./bootstrap --prefix=/usr  --no-system-libs --parallel=$(/usr/bin/getconf _NPROCESSORS_ONLN)
mkdir ../cmake-install
make DESTDIR="../cmake-install" install
cd ..
rm -r cmake-2.8.12.1
mv cmake-install/usr cmake-2.8.12.1
rm -r cmake-install

echo ""
echo ""
echo ""
echo ""
if ls cmake-2.8.12.1/bin/cmake 2>/dev/null 1>/dev/null; then
    echo "Building cmake was successful!"
else
    echo "Encountered a problem while building cmake!"
    echo "This script is not smart enough fix it :("
fi


