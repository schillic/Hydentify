
VERSION=R2014a

g++ -c  -I/usr/local/include -I/usr/local/MATLAB/$VERSION/extern/include -I/usr/local/MATLAB/$VERSION/simulink/include  -DMATLAB_MEX_FILE  -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "ppl_mexif.cpp"

g++ -c  -I/usr/local/include -I/usr/local/MATLAB/$VERSION/extern/include -I/usr/local/MATLAB/$VERSION/simulink/include  -DMATLAB_MEX_FILE  -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "ppl_wrap.cpp"

g++ -O -pthread -shared -Wl,--version-script,/usr/local/MATLAB/$VERSION/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "pplmex.mexa64"   "ppl_mexif.o"  "ppl_wrap.o"  -lgmpxx -lppl -lgmp -Wl,-rpath-link,/usr/local/MATLAB/$VERSION/bin/glnxa64 -L/usr/local/MATLAB/$VERSION/bin/glnxa64 -lmx -lmex -lmat -lm