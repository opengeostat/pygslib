# Install Python source files for pygslib subpackages

py.install_sources(
  '__init__.py',
  'version.py',
  'progress.py',
  'surpac.py',
  'utils.py',
  subdir: 'pygslib',
)

# Subpackages
subdir('gslib')
subdir('plothtml')
subdir('charttable')

# Install data directory
install_subdir(
  'data',
  install_dir: py.get_install_dir() / 'pygslib',
  install_tag: 'python-runtime',
)
