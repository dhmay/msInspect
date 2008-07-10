setlocal
set PATH=%JAVA_HOME%\bin;%PATH%
set SRC=%CD%
set DST=%TEMP%\viewer.tmp
rd /s /q %DST%
md %DST%
md %DST%\main
md %DST%\lib

REM copy source jars

REM UNDONE: shouldn't need forms_rt and swixml

pushd %DST%\lib
copy %SRC%\commons-beanutils.jar
copy %SRC%\commons-logging.jar
copy %SRC%\commons-lang.jar
copy %SRC%\commons-collections.jar
copy %SRC%\forms_rt.jar
copy %SRC%\jcommon-1.0.0.jar
copy %SRC%\jfreechart-1.0.0.jar
copy %SRC%\jai_codec.jar
copy %SRC%\jai_core.jar
copy %SRC%\Jama-1.0.1.jar
copy %SRC%\jdom.jar
copy %SRC%\log4j.jar
copy %SRC%\mlibwrapper_jai.jar
copy %SRC%\piccolo.jar
copy %SRC%\ui.jar
copy %SRC%\xercesImpl.jar
copy %SRC%\xml-apis.jar
copy %SRC%\stax-api-1.0.jar
copy %SRC%\wstx-lgpl-2.8.jar

cd %DST%\main
copy %SRC%\viewer.jar main.jar
cd %DST%

REM expand one-jar

jar -xf %SRC%\one-jar-boot.jar
copy %SRC%\one-jar.properties

REM make applicaiton jar

jar -cfm viewerApp.jar boot-manifest.mf .
move viewerApp.jar %SRC%\

popd
