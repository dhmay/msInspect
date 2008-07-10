#!/bin/sh
# set PATH=%JAVA_HOME%\bin;%PATH%
SRC=`pwd`
DST=/tmp/viewer.tmp
rm -rf $DST
mkdir $DST
mkdir $DST/main
mkdir $DST/lib

# copy source jars

# UNDONE: shouldn't need forms_rt and swixml

cd $DST/lib
cp $SRC/commons-beanutils.jar .
cp $SRC/commons-lang.jar .
cp $SRC/commons-logging.jar .
cp $SRC/commons-collections.jar .
cp $SRC/forms_rt.jar .
cp $SRC/jcommon-1.0.0.jar .
cp $SRC/jfreechart-1.0.0.jar .
cp $SRC/jai_codec.jar .
cp $SRC/jai_core.jar .
cp $SRC/Jama-1.0.1.jar .
cp $SRC/jdom.jar .
cp $SRC/log4j.jar .
cp $SRC/mlibwrapper_jai.jar .
cp $SRC/piccolo.jar .
cp $SRC/ui.jar .
cp $SRC/stax-api-1.0.jar .
cp $SRC/wstx-lgpl-2.8.jar .
cp $SRC/xercesImpl.jar .
cp $SRC/xml-apis.jar .

cd $DST/main
cp $SRC/viewer.jar main.jar
cd $DST

# expand one-jar

jar -xf $SRC/one-jar-boot.jar
cp $SRC/one-jar.properties .
# dmay adding for bugfix
rm META-INF/MANIFEST.MF

# make applicaiton jar

jar -cfm viewerApp.jar boot-manifest.mf .
mv viewerApp.jar $SRC

cd $SRC
rm -rf $DST
