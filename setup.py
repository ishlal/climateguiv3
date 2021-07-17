from cx_Freeze import setup, Executable

setup(name = 'promethee-guiv2',
      version = '0.1',
      description = 'Climate GUI',
      executables = [Executable("promethee-guiv2.py")])