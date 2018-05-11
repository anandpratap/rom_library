import os
build_examples = True

env = Environment(ENV=os.environ)
libraries = ["src", "c_api", "fortran_api"]

arma_dir = os.environ["ARMA_DIR"]
env.Append(CPPPATH=[os.path.join(arma_dir, "include")])
env.Append(LIBPATH=[os.path.join(arma_dir, "lib")])

env.Append(CXXFLAGS=["-DARMA_DONT_USE_HDF5"])
env.Append(CXXFLAGS=["-std=c++11", "-Wall", "-Wextra"])

env['CC'] = env["CXX"]
env['LINK'] = env["CXX"]

for library in libraries:
    env.SConscript(os.path.join(library, "SConstruct"), exports="env", variant_dir=os.path.join("#build", library), duplicate=0)

env.Append(LIBPATH=["#/lib"])

if build_examples:
    env.SConscript(os.path.join("examples", "SConstruct"), exports="env", variant_dir=os.path.join("#build", "examples"), duplicate=0)
