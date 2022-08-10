# example modified from https://numpy.org/devdocs/f2py/buildtools/distutils.html


from numpy.distutils.core import Extension

ext1 = Extension(name = 'add',
                 sources = ['add.f90'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'add',
          description       = "test compiling simple example via numpy.distutils",
          ext_modules = [ext1,]
          )
