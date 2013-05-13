# vim:syn=python
# encoding: ISO8859-1

APPNAME = 'divvy'
VERSION = '0.9'

def options(ctx):
  ctx.load('compiler_c')

def configure(conf):
  conf.load('compiler_c')
  conf.find_program('mpicc')
  conf.check_cfg(path='mpicc', package='', uselib_store='mpi', args='-show')
  conf.check_cfg(package='libpcre', args='--cflags --libs', uselib_store='pcre')

def build(bld):
  bld.program(source='divvy.c', target='divvy', use='mpi pcre')

