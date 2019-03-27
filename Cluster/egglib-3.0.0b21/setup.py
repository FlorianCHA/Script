import os, sys, re, glob, shutil
from distutils.core import setup, Extension
from distutils.dir_util import remove_tree
from distutils.command.clean import clean
import distutils.cmd

if '-h' in sys.argv[1:] or 'help' in sys.argv[1:]:

    sys.exit("""Welcome to EggLib's setup script

To see the manual for standard distutils commands, type: `python setup.py --help`

This page describes features specific to EggLib.

Usage:
    python setup.py build <options>
    python setup.py clean
    python setup.py uninstall

Options for the `build` command:
    --idir=<path> ........ add an include path at compilation
    --ldir=<path> ........ add a link path at link step
    --lname=<path> ....... add a library name at link step
    --cflag=<option> ..... add a compiler option
    --lflag=<option> ..... add a linker option
    --modname=<name> ..... change the name of the module
    --debug .............. add -g -O0 -DDEBUG compiler flags

    Options `--idir`, `--ldir`, `--lname`, `--cflag`, and `--lflag` can
    be repeated any number of times.

""")

# parse arguments

CPATHS = []
LPATHS = []
LNAMES = []
CFLAGS = []
LFLAGS = []
ENABLE_ABC = False
NAME = 'egglib'

argv = []
for arg in sys.argv[1:]:
    o1 = re.match('^(--idir)=(.+)$', arg)
    o2 = re.match('^(--ldir)=(.+)$', arg)
    o3 = re.match('^(--lname)=(.+)$', arg)
    o4 = re.match('^(--cflag)=(.+)$', arg)
    o5 = re.match('^(--lflag)=(.+)$', arg)
    o6 = re.match('^(--modname)=([a-zA-Z_][a-zA-Z0-9_]+)$', arg)
    o7 = re.match('^--debug$', arg)
    if o1: CPATHS.append(o1.group(2))
    elif o2: LPATHS.append(o2.group(2))
    elif o3: LNAMES.append(o3.group(2))
    elif o4: CFLAGS.append(o4.group(2))
    elif o5: LFLAGS.append(o5.group(2))
    elif o6: NAME =o6.group(2)
    elif o7: CFLAGS.extend(('-g', '-O0', '-DDEBUG'))
    else: argv.append(arg)
    if o6 and 'build' not in sys.argv[1:]:
        sys.exit('option `{0}` is available only for the `build` command'.format(res.group(1)))
sys.argv[1:] = argv

if 'install' in sys.argv[1:] or 'uninstall' in sys.argv[1:]:
    if os.path.isfile('.egglib_name'):
        NAME = open('.egglib_name').read()

if 'build' in sys.argv[1:]:
    f = open('.egglib_name', 'w')
    f.write(NAME)
    f.close()

# get the extension module files
sources = glob.glob('cppfiles/*.cpp')

# set proper options in case of linking to the GSL
if not ENABLE_ABC:
    macros = [('__woABC', '1')]
else:
    macros = []
    LNAMES.extend(['gsl', 'gslcblas'])

# make special commands
class clean_all(clean):
    def run(self):
        c = clean(self.distribution)
        c.all = True
        c.finalize_options()
        c.run()

class uninstall(distutils.cmd.Command):
    description = 'remove egglib from python distribution'
    user_options = [('name=', None, 'name of installed module')]
    def initialize_options(self):
        self.name = None
    def finalize_options(self):
        if self.name is None:
            self.name = NAME
    def run(self):
        print 'module name:', self.name
        try: mod = __import__(self.name)
        except ImportError: sys.exit('the package is not installed')
        checkfiles = os.listdir('pyfiles')
        for i in checkfiles:
            if os.path.splitext(i)[1] == '.py':
                checkfiles.append(i+'c')
        checkfiles.append('__eggwrapper.so')
        for item in os.listdir(mod.__path__[0]):
            if not item in checkfiles:
                sys.exit('this file is not expected: {0} (cancel uninstall)'.format(item))
        print 'deleting', mod.__path__[0]
        shutil.rmtree(mod.__path__[0])
        for fname in glob.iglob(mod.__path__[0] + '-*-py*.egg-info'):
            print 'removing', fname
            os.unlink(fname)

# define the extension module
binding = Extension(NAME + '._eggwrapper',
                    sources = sources,
                    include_dirs = CPATHS,
                    define_macros=macros,
                    extra_compile_args = ['-O3', '-g0'] + CFLAGS,
                    extra_link_args = LFLAGS,
                    libraries = LNAMES,
                    library_dirs = LPATHS)

# run the setup
setup(  name = NAME,
        cmdclass={'clean': clean_all, 'uninstall': uninstall},
        url = 'mycor.nancy.inra.fr/egglib/',
        author = 'Stephane De Mita, Mathieu Siol',
        author_email = 'demita@gmail.com',
        license = 'GPL v3',
        version = '3.0.0b21',                                 #EGGVERSION#
        package_dir = {NAME: 'pyfiles'},
        package_data = {NAME: ['wrappers/apps.conf']},
        packages = [NAME, NAME+'.coalesce', NAME+'.tools',
                     NAME+'.stats', NAME+'.io', NAME+'.fit',
                     NAME+'.wrappers'],
        ext_modules = [binding],
        platforms = ['cross-platform', 'Windows i686', 'MacOSX AMD64'],
        description = 'Evolutionary Genetics and Genomics Library',
        long_description = 'EggLib is a C++/Python library and program package for evolutionary genetics and genomics. Main features are sequence data management, sequence polymorphism analysis, and coalescent simulations. EggLib is a flexible Python module with a performant underlying C++ library and allows fast and intuitive development of Python programs and scripts.'
          )
