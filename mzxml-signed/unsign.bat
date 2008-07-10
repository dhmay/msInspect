@SETLOCAL
SET FILE=%1
rmdir /S /Q temp
mkdir temp
copy %FILE% temp
cd temp
jar xvf %FILE%
del %FILE%
rmdir /S /Q META-INF
jar cvf %FILE% .
cd ..
ren %FILE% %FILE%.backup
move temp\%FILE% %FILE%
rmdir /S /Q temp
