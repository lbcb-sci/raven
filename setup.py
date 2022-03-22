import os
import shutil
import subprocess

from pathlib import Path
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

from distutils import log

EXTENSION_NAME = "ravenpy"

class CMakeExtension(Extension):
  def __init__(self, name, sourcedir=""):
    Extension.__init__(self, name, sources=[])
    self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
  def build_extension(self, ext):
    install_dir = Path(self.get_ext_fullpath(ext.name)).absolute().parent
    build_dir = install_dir.joinpath('ravenpy_build')
    src_dir = Path(__file__).absolute().parent

    self.announce(f'[ravenpy::setup.py] extdir: {install_dir}', level=log.INFO)

    cmake_args = [
      f'-S {src_dir}',
      f'-B {build_dir}',

      f'-DRAVEN_BUILD_PYTHON=1',
      f'-DCMAKE_BUILD_TYPE=Release',

      '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
    ]

    cmake_args += ["-GNinja"]
    
    build_args = []
    if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
    # self.parallel is a Python 3 only way to set parallel jobs by hand
    # using -j in the build_ext call, not supported by pip or PyPA-build.
      if hasattr(self, "parallel") and self.parallel:
          # CMake 3.12+ only.
          build_args += [f"-j{self.parallel}"]

    if not os.path.exists(self.build_temp):
        os.makedirs(self.build_temp)

    fullname = self.get_ext_fullname(ext.name)
    filename = self.get_ext_filename(fullname)
  
    filepath_src = build_dir.joinpath(f'lib/{filename}') 
    filepath_dst = install_dir.joinpath(filename)
    self.announce(f'[ravenpy::setup.py] {filepath_src} -> {filepath_dst}\n', level=log.INFO)

    subprocess.check_call(
        ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)
    subprocess.check_call(
        ["cmake", "--build", build_dir] + build_args, cwd=self.build_temp)

    shutil.copy(filepath_src, filepath_dst)

if __name__ == '__main__':
  setup(
    name=EXTENSION_NAME,
    version='1.8.1',
    maintainer='Tvrtko Brekalo',
    author='Robert Vaser',
    platforms=['linux'],
    license='mit',

    ext_modules=[CMakeExtension(EXTENSION_NAME)],
    cmdclass= {
      'build_ext': CMakeBuild
    }
  )
