; msInspect.nsi
;
; This script is based on example1.nsi, but it remembers the directory, 
; has uninstall support and (optionally) installs start menu shortcuts.

; It will install msInspect into a directory that the user selects,

;--------------------------------

!define MSINSPECT_DIR "z:\cpl\msinspect"
!define MSINSPECT_DIST_DIR "${MSINSPECT_DIR}\build\dist"
!define MSINSPECT_WINDOWS_DIR "${MSINSPECT_DIR}\windows"

 !define MUI_ICON "${MSINSPECT_WINDOWS_DIR}\cpl_favicon.ico"
!define PRODUCT_STARTMENU_REGVAL "NSIS:StartMenuDir"
  Var STARTMENU_FOLDER

!include "MUI.nsh"
  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER
  !insertmacro MUI_PAGE_INSTFILES
    # These indented statements modify settings for MUI_PAGE_FINISH
    !define MUI_FINISHPAGE_NOAUTOCLOSE
    !define MUI_FINISHPAGE_SHOWREADME_CHECKED
    !define MUI_FINISHPAGE_SHOWREADME $INSTDIR\README.txt

  !insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"






; The name of the installer
Name "msInspect"

; The file to write
OutFile "msInspect_setup.exe"

; The default installation directory
InstallDir $PROGRAMFILES\msInspect

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\msInspect" "Install_Dir"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages

#Page components
#Page directory
#Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

Function .onInit
	# the plugins dir is automatically deleted when the installer exits
	InitPluginsDir
	File /oname=$PLUGINSDIR\splash.bmp "${MSINSPECT_WINDOWS_DIR}\splash.bmp"

	splash::show 1000 $PLUGINSDIR\splash

	Pop $0 ; $0 has '1' if the user closed the splash screen early,
			; '0' if everything closed normally, and '-1' if some error occurred.
FunctionEnd


;--------------------------------

; The stuff to install
Section "msInspect (required)"

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  ; Put file there
  ;  !insertmacro ZIPDLL_EXTRACT "$INSTDIR\viewer-dist.zip" $INSTDIR "<ALL>"
  File 	${MSINSPECT_DIST_DIR}\*.*
  File  ${MSINSPECT_WINDOWS_DIR}\README.txt
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\msInspect "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\msInspect" "DisplayName" "msInspect"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\msInspect" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\msInspect" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\msInspect" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
SectionEnd


; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    CreateShortCut "$DESKTOP\msInspect.lnk" "$INSTDIR\msInspect.exe"
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\msInspect.lnk" "$INSTDIR\msInspect.exe"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\README.lnk" "notepad.exe" "$INSTDIR\README.txt"
  !insertmacro MUI_STARTMENU_WRITE_END
SectionEnd





;--------------------------------

; Uninstaller

Section "Uninstall"
  
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\msInspect"
  DeleteRegKey HKLM SOFTWARE\msInspect

  ; Remove files and uninstaller
  Delete $INSTDIR\*.*

  ; Remove shortcuts, if any
  Delete "$SMPROGRAMS\msInspect\*.*"

  ; Remove desktop shortcut
  Delete "$DESKTOP\msInspect.lnk"

  ; Remove directories used
  RMDir "$SMPROGRAMS\msInspect"
  RMDir "$INSTDIR"

SectionEnd
