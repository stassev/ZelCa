#!/bin/bash

if ls install.env 2>/dev/null 1>/dev/null; then
    echo "Checking install.env ..."
    if grep EIGENINCLUDES install.env 1>/dev/null 2>/dev/null; then 
        if (ls $(grep EIGENINCLUDES install.env | tail -n 1  | awk -F"=" '{print $2}')/Eigen/LU  2>/dev/null 1>/dev/null) ||  (ls DEPENDENCIES/$(grep EIGENINCLUDES install.env | tail -n 1  | awk -F"=" '{print $2}')/Eigen/LU  2>/dev/null 1>/dev/null); then
            echo "Eigen found."
        else  
            echo "Missing Eigen!"
            echo "Run this again and tell me its location:"
            echo "./build.sh initial-zelca"
            echo "or if you want to build it directly then run:"
            echo "./build.sh eigen"
            exit 1
        fi
    else
        echo "Not sure where Eigen is installed."
        echo "Run this again and tell me:"
        echo "./build.sh initial-zelca"
        echo "or if you want to build it directly then run:"
        echo "./build.sh eigen"
        exit 1
    fi
    
    if grep CUBAINCLUDES install.env 1>/dev/null 2>/dev/null; then 
        if (ls $(grep CUBAINCLUDES install.env | tail -n 1  | awk -F"=" '{print $2}')/cuba.h  2>/dev/null 1>/dev/null) ||  (ls DEPENDENCIES/$(grep CUBAINCLUDES install.env | tail -n 1  | awk -F"=" '{print $2}')/cuba.h  2>/dev/null 1>/dev/null); then
            echo "cuba.h found."
        else  
            echo "Missing cuba.h!"
            echo "Run this again and tell me its location:"
            echo "./build.sh initial-zelca"
            echo "or if you want to build it directly then run:"
            echo "./build.sh cuba"
            exit 1
        fi
    else
        echo "Not sure where cuba.h is installed."
        echo "Run this again and tell me:"
        echo "./build.sh initial-zelca"
        echo "or if you want to build it directly then run:"
        echo "./build.sh cuba"
        exit 1
    fi
    
    if grep CUBALIBS install.env 1>/dev/null 2>/dev/null; then 
        if (ls $(grep CUBALIBS install.env | tail -n 1  | awk -F"=" '{print $2}')/libcuba.a  2>/dev/null 1>/dev/null) ||  (ls DEPENDENCIES/$(grep CUBALIBS install.env | tail -n 1  | awk -F"=" '{print $2}')/libcuba.a  2>/dev/null 1>/dev/null); then
            echo "libcuba.a found."
        else  
            echo "Missing libcuba.a!"
            echo "Run this again and tell me its location:"
            echo "./build.sh initial-zelca"
            echo "or if you want to build it directly then run:"
            echo "./build.sh cuba"
            exit 1
        fi
    else
        echo "Not sure where libcuba.a is installed."
        echo "Run this again and tell me:"
        echo "./build.sh initial-zelca"
        echo "or if you want to build it directly then run:"
        echo "./build.sh cuba"
        exit 1
    fi
else
    echo "install.env is missing!"
    echo "You first need to run:"
    echo "./build.sh initial-zelca"
    exit 1
fi
echo ""
echo "******************************** BEGIN ZelCa BUILD ********************************"
echo "******************************** BEGIN ZelCa BUILD ********************************"
echo ""

export INCLUDES=" -I$(grep EIGENINCLUDES install.env | tail -n 1  | awk -F"=" '{print $2}') "
INCLUDES+=" -I$(grep CUBAINCLUDES install.env |  tail -n 1 | awk -F"=" '{print $2}') "
export LIBS=" -L$(grep CUBALIBS install.env  | tail -n 1 | awk -F"=" '{print $2}') "

cd ZelCa-source
rm ../executable/ZelCa
make clean
make
make clean
mv ZelCa ../executable/
echo ""
echo ""
echo ""
echo ""
if ls ../executable/ZelCa 2>/dev/null 1>/dev/null; then
    echo "Zelca is successfully built!"
    echo "The executable is located at: ./executable/Zelca"
else
    echo "Something broke!"
    echo "This script is not smart enough to fix it :("
fi
cd ..


