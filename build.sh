#!/bin/bash

case "$1" in 
    cmake)
        ./scripts/cmake-build.sh
        ;;
        
    cuba)
        ./scripts/cuba-build.sh
        ;;
        
    eigen)
        ./scripts/eigen-build.sh
        ;;
    
    zelca)
        ./scripts/zelca-build.sh
        ;;
    
    CLEAN-ALL)
        rm -r DEPENDENCIES 2>/dev/null
        rm executable/ZelCa 2>/dev/null
        rm ZelCa-source/*.o  2>/dev/null
        rm ZelCa-source/ZelCa  2>/dev/null
        rm install.env 2>/dev/null
        ;;
    
    initial-zelca)
        echo ""
        echo ""
        echo ""
        echo "Welcome to the initial setup of ZelCa!"
        echo "ZelCa needs some (yet) not widely used libraries: eigen (v3.2) and Cuba (v3.2). Eigen (http://eigen.tuxfamily.org) is a C++ template library for linear algebra. Cuba (http://www.feynarts.de/cuba/) is a library for multidimensional numerical integration. " 
        echo "This setup will guide you through building those libraries (if missing) on your machine before attempting to build ZelCa."
        echo "Note that Eigen needs the cmake build automation software. Cmake seems to be missing from some clusters. Therefore, if you need to build Eigen, this setup will attempt to build cmake as well (if missing)."
        echo "Note that only the contents of this directory will be modified. All dependencies will be built in a subdirectory called DEPENDENCIES."
        echo "Let's begin!"
        echo ""
        echo ""
        echo ""
        echo "Do you wish to (re)build Cuba-3.2? If not, I will use your system-wide Cuba if present." 
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) ./build.sh cuba; break;;
                No ) echo "Will try to use system-wide installed Cuba-3.2."
                     if ls /usr/lib/libcuba.a 2>/dev/null 1>/dev/null; then
                        echo "CUBALIBS=/usr/lib/" >> install.env
                     else 
                        read -p "Please enter the path containing libcuba.a (e.g. /usr/lib/): " cubalibDir
                        echo "CUBALIBS=${cubalibDir}/" >> install.env
                     fi   
                     if ls /usr/include/cuba.h 2>/dev/null 1>/dev/null; then
                        echo "CUBAINCLUDES=/usr/include/" >> install.env
                     else 
                        read -p "Please enter the path containing cuba.h (e.g. /usr/include/): " cubahDir
                        echo "CUBAINCLUDES=${cubahDir}/" >> install.env
                     fi   
                     break;;
            esac
        done
        echo ""
        echo ""
        echo ""
        echo "Do you wish to (re)build Eigen-3.2? If not, I will use your system-wide Eigen if present."
        select yn in "Yes" "No"; do
            case $yn in
                Yes )  echo ""
                       echo ""
                       echo ""
                       echo "First we need cmake."
                       echo "Do you wish to (re)build cmake? If not, I will use your system-wide cmake if present." 
                       select yn in "Yes" "No"; do
                            case $yn in
                                Yes ) ./build.sh cmake; break;;
                                No ) echo "Will try to use system-wide installed cmake."; break;;
                            esac
                       done
                       ./build.sh eigen
                       break;; 
                No ) echo "Will try to use system-wide installed Eigen-3.2."
                     if ls /usr/include/eigen3/Eigen/LU 2>/dev/null 1>/dev/null; then
                        echo "EIGENINCLUDES=/usr/include/eigen3/" >> install.env
                     else 
                        echo "Please enter the path to your Eigen installation "
                        read -p "(e.g. enter /usr/include/eigen3/ if file LU is in /usr/include/eigen3/Eigen/): " eigenDir
                        echo "EIGENINCLUDES=${eigenDir}/" >> install.env
                     fi   
                     break;;
            esac
        done
        ./build.sh zelca
        ;;
        
    *)
        echo "Usage:"
        echo "-- Initial build:"
        echo "      $0 initial-zelca"
        echo "   The above command should guide you in building all missing dependencies."
        echo "   If you want to build any of the depencies at a later time then choose:"
        echo "      $0 cuba"
        echo "      $0 eigen"
        echo "      $0 cmake"
        echo "-- Subsequent builds (assumes all dependencies are met):"
        echo "      $0 zelca"
        echo "-- To return to prestine state:"
        echo "      $0 CLEAN-ALL"
        exit 1
        
esac


