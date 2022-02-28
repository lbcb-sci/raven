from distutils.errors import DistutilsExecError
from typing import List, Optional

from collections.abc import Mapping

import pathlib
import shutil
import os

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

from distutils import log

EXTENSION_NAME = "RavenPy"

class CMakeExtension(Extension):
  ''' Addapter class for setuptools Extension.
      Used for naming; assumes source files are handled by CMakeBuildExt'''

  def __init__(self, name: str):
    super().__init__(name=name, sources=[])

class CMakeBuildExt(build_ext):
  ''' Executes cmake build procedure '''

  def run(self) -> None:
    for ext in self.extensions:
      if ext.name == EXTENSION_NAME:
        self.__cmake_build(ext)

  @staticmethod
  def __find_cmake() -> Optional[str]:
    return shutil.which('cmake')

  def __cmake_build(self, ext):
    cmake_path_str: Optional[str] = self.__find_cmake()
    if cmake_path_str is None:
      raise DistutilsExecError('[RavenPy::setup.py] failed to locate cmake')

    base_build_dir = pathlib.Path(self.build_temp).absolute()

    ravenpy_build_path: pathlib.Path = base_build_dir.joinpath('ravenpy_build')
    ravenpy_install_path: pathlib.Path = base_build_dir.joinpath('ravenpy_install')

    os.makedirs(base_build_dir, exist_ok=True)
    self.announce(f'[RavenPy::setup.py] build directory {ravenpy_build_path}', level=log.INFO)
    self.announce(f'[RavenPy::setup.py] install directory {ravenpy_install_path}', level=log.INFO)


    build_type: str = 'Debug' if self.debug else 'Release'

    install_libdir: str = 'lib'
    install_includedir = 'include'

    cmake_config: List[str] = [
      '-H./',
      f'-B{ravenpy_build_path}',
      f'-DCMAKE_BUILD_TYPE={build_type}',

      f'-DRAVEN_BUILD_PYTHON=1',

      f'-DCMAKE_INSTALL_PREFIX={ravenpy_install_path}',
      f'-DCMAKE_INSTALL_INCLUDEDIR={install_includedir}',
      f'-DCMAKE_INSTALL_LIBDIR={install_libdir}'
    ]

    try:
      self.announce(f'[RavenPy::setup.py] generating cmake files', level=log.INFO)
      self.spawn([cmake_path_str, *cmake_config])
    except DistutilsExecError as dee:
      self.announce(f'[RavenPy::setup.py] failed to configure cmake', level=log.FATAL)
      raise dee

    cmake_build_args: List[str] = [
        '--build', str(ravenpy_build_path),
        '--config', build_type,
        '--target', 'install'
      ]

    try:
      self.announce(f'[RavenPy::setup.py]')
      self.spawn([cmake_path_str, *cmake_build_args])
    except DistutilsExecError as dee:
      self.announce(f'[RavenPy::setup.py] failed to build binaries', level=log.FATAL) 
      raise dee

if __name__ == '__main__':
  setup(
    name=EXTENSION_NAME,
    version='1.8.0',
    maintainer='Tvrtko Brekalo',
    author='Robert Vaser',
    platforms=['linux'],
    license='mit',

    ext_modules=[CMakeExtension(EXTENSION_NAME)],
    cmdclass= {
      'build_ext': CMakeBuildExt
    }
  )
