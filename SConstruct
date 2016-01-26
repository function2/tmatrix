# Build instructions.
# Run 'scons' to build.

import os,os.path,re,subprocess

# Hide build messages. Comment this for verbose output.
# env.SetDefault(
#     CCCOMSTR='Compiling $TARGET',
#     CXXCOMSTR='Compiling $TARGET',
#     LINKCOMSTR='Linking $TARGET',
# )

BUILD_DIR = 'build'
SRC_DIR = 'src'
OBJ_DIR = '%s/obj' % (BUILD_DIR)

# Select compiler
# CXX = 'g++' , CXX = 'clang++'
env = Environment(CXX = 'eg++')

# Put scons database in build directory.
env.SConsignFile("%s/scons/scons-signatures" % BUILD_DIR)

# Flags.
include_flags = ' -I./include -I/usr/local/include '
BUILD_FLAGS = '-O9 -Wall' + include_flags
link_flags = '-s -L/usr/local/lib '
GMPXX_FLAGS = R' -DMATRIX_GMPXX -DFLOAT_ZERO_EPS=\\\"1e-16\\\" -lgmpxx -lgmp'
env.Append(LINKFLAGS = link_flags)

# comment this out for older compilers.
BUILD_FLAGS += ' -std=c++11'

# File names.
input_names = [
    'LDU',
    'check_eigenval',
    'cofactor',
    'determinant',
    'exponent',
    'inverse',
    'is_normal',
    'mult',
    'pow',
    'proj',
    'qr',
    'rref',
    'subspaces',
    'x_r',
    'positive_definite',
    'sandbox',
]

# Files to build with rational numbers only.
rat_names = [
    'chem_eq_balance',
]

# Files to build with floating point only.
float_names = [
    'eigen2',
]

def build_it(target,source,object,flags):
    srcfile = '%s/' % (SRC_DIR) + source + '.cpp'
    objfile = '%s/' % (OBJ_DIR) + object + '.o'
    env.StaticObject(objfile,srcfile, parse_flags = flags)
    return env.Program('%s/' % (BUILD_DIR) + target, objfile, parse_flags = flags)

# Build binaries using rational numbers.
rational_flags = BUILD_FLAGS + GMPXX_FLAGS + ' -DMATRIX_RATIONAL'
for g in input_names + rat_names:
    build_it(g,g,g + "_rational",rational_flags)

# Build binaries using floating points.
float_flags = BUILD_FLAGS + GMPXX_FLAGS + ' -DMATRIX_FLOAT -DMATRIX_FLOAT_PREC=1024'
for g in input_names + float_names:
    build_it(g + '_float',g,g + '_float', float_flags)
