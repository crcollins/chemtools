# -*- mode: python -*-
a = Analysis([os.path.join(HOMEPATH,'support/_mountzlib.py'), os.path.join(HOMEPATH,'support/useUnicode.py'), 'fileparser.py'],
             pathex=['/home/chris/Desktop/chemtools/fileparser'])
pyz = PYZ(a.pure)
exe = EXE( pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name=os.path.join('dist', 'fileparser'),
          debug=False,
          strip=False,
          upx=True,
          console=1 )
