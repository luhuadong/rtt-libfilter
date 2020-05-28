from building import *
Import('rtconfig')

src   = []
cwd   = GetCurrentDir()

# add libfilter src files.
if GetDepend('PKG_USING_LIBFILTER'):
    src += Glob('src/xxx.c')

if GetDepend('PKG_USING_LIBFILTER_SAMPLE'):
    src += Glob('examples/xxx_sample.c')

# add libfilter include path.
path  = [cwd + '/inc']

# add src and include to group.
group = DefineGroup('libfilter', src, depend = ['PKG_USING_LIBFILTER'], CPPPATH = path)

Return('group')
