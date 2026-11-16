Option Explicit

Dim fso, shell, baseDir, serverScript, localAppData, pyExe, candidates, cand
Set fso = CreateObject("Scripting.FileSystemObject")
Set shell = CreateObject("WScript.Shell")

baseDir = fso.GetParentFolderName(WScript.ScriptFullName)
serverScript = fso.BuildPath(baseDir, "server_silent.py")
If Not fso.FileExists(serverScript) Then WScript.Quit 1

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
shell.Run """" & pyExe & """ """ & serverScript & """ --stop", 0, True
If Err.Number <> 0 Then
  Err.Clear
  shell.Run "python """ & serverScript & """ --stop", 0, True
End If
On Error GoTo 0

