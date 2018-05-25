import os
AddOption("--build-examples", action="store_true", dest="build_examples", default=False)
AddOption("--intel-mpi", action="store_true", dest="intel_mpi", default=False)
AddOption("--intel", action="store_true", dest="intel", default=False)
AddOption("--dbg", action="store_true", dest="dbg", default=False)

env = Environment(ENV=os.environ)
libraries = ["src", "c_api", "fortran_api", "python_api"]

arma_dir = os.environ["ARMA_DIR"]
env.Append(CPPPATH=[os.path.join(arma_dir, "include")])

if os.path.isdir(os.path.join(arma_dir, "lib")):
    env.Append(LIBPATH=[os.path.join(arma_dir, "lib")])
elif os.path.isdir(os.path.join(arma_dir, "lib64")):
    env.Append(LIBPATH=[os.path.join(arma_dir, "lib64")])
   
env.Append(CXXFLAGS=["-DARMA_DONT_USE_HDF5"])
env.Append(CCFLAGS=["-DARMA_DONT_USE_HDF5"])
env.Append(CXXFLAGS=["-std=c++11", "-Wall", "-Wextra"])

if GetOption("intel_mpi"):
    env['CXX'] = "mpicxx"
    env['F90'] = "mpiifort"
    env.Append(F90FLAGS=["-cpp"])
elif GetOption("intel"):
    env['CXX'] = "icpc"
    env['F90'] = "ifort"
    env.Append(F90FLAGS=["-cpp"])

if GetOption("dbg"):
    env.Append(F90FLAGS=["-g"])
    env.Append(CCFLAGS=["-g"])
    env.Append(CXXFLAGS=["-g"])

env['CC'] = env["CXX"]
env['LINK'] = env["F90"]


for library in libraries:
    env.SConscript(os.path.join(library, "SConstruct"), exports="env", variant_dir=os.path.join("#build", library), duplicate=0)

env.Append(LIBPATH=["#/lib"])

if GetOption("build_examples"):
    env.SConscript(os.path.join("examples", "SConstruct"), exports="env", variant_dir=os.path.join("#build", "examples"), duplicate=0)

combine_libs = Builder(action = "ar -rcT $TARGET $SOURCES")
env.Append(BUILDERS = {'combine_libs' : combine_libs})
env.combine_libs("lib/librom.a", Glob("lib/libmain_*.a"))
                       
