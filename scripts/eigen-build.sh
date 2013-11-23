#!/bin/bash

mkdir DEPENDENCIES 2>/dev/null
cd DEPENDENCIES

if ls cmake-2.8.12.1/bin/cmake 2>/dev/null 1>/dev/null; then
    if hash cmake 2>/dev/null; then
        echo "You have a local cmake in DEPENDENCIES"
        echo "as well as a system-wide installed cmake."
        echo "Select which one to use to build eigen:"
        select cm in "local" "system-wide"; do
             case $cm in
                 "local" )       CMAKE=../cmake-2.8.12.1/bin/cmake
                                 break;;
                 "system-wide" ) CMAKE=cmake
                                 break;;
             esac
        done
    else
        CMAKE=../cmake-2.8.12.1/bin/cmake
    fi
else
    if hash cmake 2>/dev/null; then
        CMAKE=cmake
    else
        echo "Encountered a problem: Missing CMAKE!"
        echo "Build cmake first by invoking:"
        echo "      ./build.sh cmake"
        exit 1
    fi
fi

echo "******************************** BEGIN EIGEN BUILD ********************************"
echo "******************************** BEGIN EIGEN BUILD ********************************"
echo "******************************** BEGIN EIGEN BUILD ********************************"
echo "******************************** BEGIN EIGEN BUILD ********************************"
echo "******************************** BEGIN EIGEN BUILD ********************************"
rm -r eigen-3.2.0* 2> /dev/null
wget https://bitbucket.org/eigen/eigen/get/3.2.0.tar.bz2 --output-document=eigen-3.2.0.tar.bz
mkdir eigen-3.2.0
tar -xvjf eigen-3.2.0.tar.bz -C eigen-3.2.0 --strip-components=1
mkdir eigen-build
cd eigen-build
${CMAKE}  ../eigen-3.2.0  -DCMAKE_BUILD_TYPE=Release  -DCMAKE_INSTALL_PREFIX=/usr
mkdir ../eigen-install
make DESTDIR="../eigen-install" install
cd ..
rm -r eigen-build
rm -r eigen-3.2.0
mv eigen-install/usr eigen-3.2.0
rm -r eigen-install

echo ""
echo ""
echo ""
echo ""
if ls eigen-3.2.0/include/eigen3/Eigen/LU 2>/dev/null 1>/dev/null; then
    echo "Building eigen was successful!"
    echo "EIGENINCLUDES=../DEPENDENCIES/eigen-3.2.0/include/eigen3/" >> ../install.env
else
    echo "Encountered a problem while building eigen!"
    echo "This script is not smart enough fix it :("
fi
