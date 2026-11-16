Option Explicit

Dim fso, shell, baseDir, distIndex, serverScript, localAppData, pyExe, candidates, cand
Set fso = CreateObject("Scripting.FileSystemObject")
Set shell = CreateObject("WScript.Shell")

baseDir = fso.GetParentFolderName(WScript.ScriptFullName)
distIndex = fso.BuildPath(baseDir, "dist\index.html")
serverScript = fso.BuildPath(baseDir, "server_silent.py")

If Not fso.FileExists(distIndex) Then
  MsgBox "dist\index.html 파일이 없습니다. 먼저 'npm run build'를 실행해주세요.", vbExclamation, "YeastWGS Studio"
  WScript.Quit 1
End If

If Not fso.FileExists(serverScript) Then
  MsgBox "server_silent.py 파일을 찾을 수 없습니다.", vbCritical, "YeastWGS Studio"
  WScript.Quit 1
End If

localAppData = shell.ExpandEnvironmentStrings("%LOCALAPPDATA%")
candidates = Array( _
  localAppData & "\Programs\Python\Python313\pythonw.exe", _
  localAppData & "\Programs\Python\Python312\pythonw.exe", _
  localAppData & "\Programs\Python\Python311\pythonw.exe", _
  localAppData & "\Programs\Python\Python310\pythonw.exe" _
)

pyExe = ""
For Each cand In candidates
  If fso.FileExists(cand) Then
    pyExe = cand
    Exit For
  End If
Next

If pyExe = "" Then
  pyExe = "pythonw"
End If

On Error Resume Next
shell.Run """" & pyExe & """ """ & serverScript & """", 0, False
If Err.Number <> 0 Then
  Err.Clear
  shell.Run "python """ & serverScript & """", 0, False
End If
On Error GoTo 0

