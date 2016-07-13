from distutils.core import setup, Extension

extension_mod = Extension("_UnsymmetricPairAlignment", sources=["UnsymmetricPairAlignment_wrap.c", "UnsymmetricPairAlignment.c"], swig_opts=['-modern', '-I../include'])

setup(name = "UnsymmetricPairAlignment", ext_modules=[extension_mod], py_modules=['UnsymmetricPairAlignment'])
