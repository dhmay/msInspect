setlocal
set DROP=%1
if "%DROP%"=="" set DROP=\msInspect
echo %DROP%

set LIBS=%CD%\..\..\external\lib\client
set XCOPY=xcopy /D

mkdir %DROP%
mkdir %DROP%\src
mkdir %DROP%\lib
mkdir %DROP%\classes

copy build-drop.xml %DROP%\build.xml
%XCOPY% log4j.xml %DROP%\src
%XCOPY% /s org %DROP%\src\org\
%XCOPY% /s com %DROP%\src\com\
%XCOPY% /s modwt %DROP%\src\modwt\

pushd %DROP%\lib
%XCOPY% %LIBS%\commons-beanutils.jar
%XCOPY% %LIBS%\commons-collections.jar
%XCOPY% %LIBS%\commons-lang.jar
%XCOPY% %LIBS%\commons-logging.jar
%XCOPY% %LIBS%\forms_rt.jar
%XCOPY% %LIBS%\j3dcore.jar
%XCOPY% %LIBS%\j3dutils.jar
%XCOPY% %LIBS%\jai_codec.jar
%XCOPY% %LIBS%\jai_core.jar
%XCOPY% %LIBS%\Jama-1.0.1.jar
%XCOPY% %LIBS%\jcommon-1.0.16.jar
%XCOPY% %LIBS%\jdom.jar
%XCOPY% %LIBS%\jfreechart-1.0.13.jar
%XCOPY% %LIBS%\junit.jar
%XCOPY% %LIBS%\log4j.jar
%XCOPY% %LIBS%\mlibwrapper_jai.jar
%XCOPY% %LIBS%\piccolo.jar
%XCOPY% %LIBS%\ui.jar
%XCOPY% %LIBS%\vecmath.jar
%XCOPY% %LIBS%\..\web\xercesImpl.jar
%XCOPY% %LIBS%\..\web\xml-apis.jar
popd

