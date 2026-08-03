project('pygslib', 'c', 'fortran', 'cython',
  version: '0.0.0.7.1.1',
  meson_version: '>=1.1.0',
  default_options: [
    'warning_level=1',
    'buildtype=release',
  ],
)

# ---------------------------------------------------------------------------
# Python
# ---------------------------------------------------------------------------
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# ---------------------------------------------------------------------------
# NumPy (>=2.0)
# ---------------------------------------------------------------------------
np_dep = dependency('numpy')

incdir_numpy = run_command(py,
  ['-c', 'import os; import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import os; import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true,
).stdout().strip()

np_compile_args = ['-I' + incdir_numpy, '-I' + incdir_f2py]

fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: np_compile_args,
  dependencies: py_dep,
  gnu_symbol_visibility: 'hidden',
)

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  compile_args: np_compile_args,
)

# ---------------------------------------------------------------------------
# Cython settings
# ---------------------------------------------------------------------------
cython_args = ['-3', '--fast-fail',
               '--output-file', '@OUTPUT@',
               '--include-dir', meson.project_source_root(),
               '@INPUT@']

cython_gen = generator(find_program('cython'),
  arguments: cython_args,
  output: '@BASENAME@.c',
)

cython_c_args = ['-DCYTHON_CCOMPLEX=0']

# ============================================================================
# FORTRAN f2py EXTENSIONS
# ============================================================================

# __rotscale — only: rotscale
f2py_rotscale = custom_target('f2py_rotscale',
  input: 'for_code/rotscale.f90',
  output: ['__rotscalemodule.c', '__rotscale-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__rotscale', '--lower',
            'only:', 'rotscale', ':'],
)
py.extension_module('__rotscale',
  ['for_code/rotscale.f90', f2py_rotscale, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __block_covariance — only: block_covariance (includes setrot, cova3, sqdist)
f2py_block_covariance = custom_target('f2py_block_covariance',
  input: 'for_code/block_covariance.f90',
  output: ['__block_covariancemodule.c', '__block_covariance-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__block_covariance', '--lower',
            'only:', 'block_covariance', ':'],
)
py.extension_module('__block_covariance',
  ['for_code/block_covariance.f90', f2py_block_covariance, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __read_gslib — only: read_header read_ndata read_data
f2py_read_gslib = custom_target('f2py_read_gslib',
  input: 'for_code/read_gslib.f90',
  output: ['__read_gslibmodule.c', '__read_gslib-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__read_gslib', '--lower',
            'only:', 'read_header', 'read_ndata', 'read_data', ':'],
)
py.extension_module('__read_gslib',
  ['for_code/read_gslib.f90', f2py_read_gslib, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __addcoord — only: addcoord
f2py_addcoord = custom_target('f2py_addcoord',
  input: 'for_code/addcoord.f90',
  output: ['__addcoordmodule.c', '__addcoord-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__addcoord', '--lower',
            'only:', 'addcoord', ':'],
)
py.extension_module('__addcoord',
  ['for_code/addcoord.f90', f2py_addcoord, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __general — only: setrot cova3 sqdist (includes setrot.f90, cova3.f90, sqdist.f90)
f2py_general = custom_target('f2py_general',
  input: 'for_code/general.f90',
  output: ['__generalmodule.c', '__general-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__general', '--lower',
            'only:', 'setrot', 'dsetrot', 'cova3', 'dcova3', 'sqdist', 'dsqdist', ':'],
)
py.extension_module('__general',
  ['for_code/general.f90', f2py_general, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __plot — only: histplt probplt qpplt (includes locate, sortem)
f2py_plot = custom_target('f2py_plot',
  input: 'for_code/plot.f90',
  output: ['__plotmodule.c', '__plot-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__plot', '--lower',
            'only:', 'histplt', 'probplt', 'qpplt', ':'],
)
py.extension_module('__plot',
  ['for_code/plot.f90', f2py_plot, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __declus — only: declus
f2py_declus = custom_target('f2py_declus',
  input: 'for_code/declus.f90',
  output: ['__declusmodule.c', '__declus-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__declus', '--lower',
            'only:', 'declus', ':'],
)
py.extension_module('__declus',
  ['for_code/declus.f90', f2py_declus, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __variograms — only: writeout gamv gamma writeout_gam gamv3d
f2py_variograms = custom_target('f2py_variograms',
  input: 'for_code/variograms.f90',
  output: ['__variogramsmodule.c', '__variograms-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__variograms', '--lower',
            'only:', 'writeout', 'gamv', 'gamma', 'writeout_gam', 'gamv3d', ':'],
)
py.extension_module('__variograms',
  ['for_code/variograms.f90', f2py_variograms, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __bigaus — only: bigaus simpsonf simpsong (includes setrot, cova3, sqdist, gauinv)
f2py_bigaus = custom_target('f2py_bigaus',
  input: 'for_code/bigaus.f90',
  output: ['__bigausmodule.c', '__bigaus-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__bigaus', '--lower',
            'only:', 'bigaus', 'simpsonf', 'simpsong', ':'],
)
py.extension_module('__bigaus',
  ['for_code/bigaus.f90', f2py_bigaus, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __bicalib — only: bicalib
f2py_bicalib = custom_target('f2py_bicalib',
  input: 'for_code/bicalib.f90',
  output: ['__bicalibmodule.c', '__bicalib-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__bicalib', '--lower',
            'only:', 'bicalib', ':'],
)
py.extension_module('__bicalib',
  ['for_code/bicalib.f90', f2py_bicalib, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __trans — only: trans (includes locate, sortem, acorni, powint, getz)
f2py_trans = custom_target('f2py_trans',
  input: 'for_code/trans.f90',
  output: ['__transmodule.c', '__trans-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__trans', '--lower',
            'only:', 'trans', ':'],
)
py.extension_module('__trans',
  ['for_code/trans.f90', f2py_trans, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __draw — only: draw (includes locate, acorni)
f2py_draw = custom_target('f2py_draw',
  input: 'for_code/draw.f90',
  output: ['__drawmodule.c', '__draw-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__draw', '--lower',
            'only:', 'draw', ':'],
)
py.extension_module('__draw',
  ['for_code/draw.f90', f2py_draw, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __dm2csv — only: dm2csv_ep dm2csv_sp
f2py_dm2csv = custom_target('f2py_dm2csv',
  input: 'for_code/dm2csv.f90',
  output: ['__dm2csvmodule.c', '__dm2csv-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@', '-m', '__dm2csv', '--lower',
            'only:', 'dm2csv_ep', 'dm2csv_sp', ':'],
)
py.extension_module('__dm2csv',
  ['for_code/dm2csv.f90', f2py_dm2csv, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# __dist_transf — only: backtr anatbl nscore ns_ttable
f2py_dist_transf = custom_target('f2py_dist_transf',
  input: 'for_code/dist_transf.f90',
  output: ['__dist_transfmodule.c', '__dist_transf-f2pywrappers.f'],
  command: [find_program('f2py'), '@INPUT@',
            '-m', '__dist_transf', '--lower',
            'only:', 'backtr', 'anatbl', 'nscore', 'ns_ttable', ':'],
)
py.extension_module('__dist_transf',
  ['for_code/dist_transf.f90', f2py_dist_transf, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# ============================================================================
# FORTRAN f2py EXTENSIONS — complex multi-file
# ============================================================================

# kt3d — 10 source files, only: pykt3d set_unest
kt3d_sources = files(
  'for_code/kt3d/kt3d.f90',
  'for_code/kt3d/gslib/setrot.f90',
  'for_code/kt3d/gslib/getindx.f90',
  'for_code/kt3d/gslib/picksupr.f90',
  'for_code/kt3d/gslib/setsupr.f90',
  'for_code/kt3d/gslib/sqdist.f90',
  'for_code/kt3d/gslib/cova3.f90',
  'for_code/kt3d/gslib/ktsol.f90',
  'for_code/kt3d/gslib/sortem.f90',
  'for_code/kt3d/gslib/srchsupr.f90',
)

f2py_kt3d = custom_target('f2py_kt3d',
  input: kt3d_sources,
  output: ['__gslib__kt3dmodule.c', '__gslib__kt3d-f2pywrappers.f', '__gslib__kt3d-f2pywrappers2.f90'],
  command: [find_program('f2py')] + kt3d_sources + [
    '-m', '__gslib__kt3d', '--lower',
    'only:', 'pykt3d', 'set_unest', ':',
  ],
)
py.extension_module('__gslib__kt3d',
  kt3d_sources + [f2py_kt3d, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# postik — 5 source files, only: postik set_unest get_unest
postik_sources = files(
  'for_code/postik/postik.f90',
  'for_code/postik/gslib/beyond.f90',
  'for_code/postik/gslib/locate.f90',
  'for_code/postik/gslib/powint.f90',
  'for_code/postik/gslib/sortem.f90',
)

f2py_postik = custom_target('f2py_postik',
  input: postik_sources,
  output: ['__gslib__postikmodule.c', '__gslib__postik-f2pywrappers.f', '__gslib__postik-f2pywrappers2.f90'],
  command: [find_program('f2py')] + postik_sources + [
    '-m', '__gslib__postik', '--lower',
    'only:', 'postik', 'set_unest', 'get_unest', ':',
  ],
)
py.extension_module('__gslib__postik',
  postik_sources + [f2py_postik, fortranobject_c],
  dependencies: [np_dep, fortranobject_dep],
  install: true, subdir: 'pygslib/gslib',
)

# ============================================================================
# CYTHON EXTENSIONS
# ============================================================================

cython_modules = {
  'drillhole':  'cython_code/drillhole.pyx',
  'blockmodel': 'cython_code/blockmodel.pyx',
  'vtktools':   'cython_code/vtktools.pyx',
  'nonlinear':  'cython_code/nonlinear.pyx',
  'sandbox':    'cython_code/sandbox.pyx',
}

foreach name, src : cython_modules
  py.extension_module(name,
    [cython_gen.process(src)],
    c_args: cython_c_args,
    include_directories: include_directories('pygslib'),
    dependencies: [np_dep, py_dep],
    install: true,
    subdir: 'pygslib',
  )
endforeach

# ============================================================================
# PURE PYTHON SOURCES
# ============================================================================

python_sources = [
  'pygslib/__init__.py',
  'pygslib/version.py',
  'pygslib/progress.py',
  'pygslib/surpac.py',
  'pygslib/utils.py',
]

py.install_sources(
  python_sources,
  subdir: 'pygslib',
)

subdir('pygslib')
